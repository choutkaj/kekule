"""HTTP access to upstream databases, used only by maintainers' selection.

Benchmark runs never call this module: they read the hash-locked snapshot.
"""
from __future__ import annotations

import hashlib
import json
import random
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

USER_AGENT = "kekule-bench-selection (+https://github.com/choutkaj/kekule)"
RETRY_STATUS = (429, 500, 502, 503, 504)


class NotFound(Exception):
    """The upstream resource does not exist (HTTP 404)."""


class Client:
    """Polite HTTP client: paced requests, retries with backoff, optional response cache.

    The cache only lets an interrupted selection resume; the published
    snapshot, not this cache, is what pins a dataset.
    """

    def __init__(self, interval: float = 0.25, attempts: int = 6, cache: Path | None = None):
        self.interval = interval
        self.attempts = attempts
        self.cache = cache
        self._lock = threading.Lock()
        self._next = 0.0

    def _pace(self, delay: float = 0.0) -> None:
        with self._lock:
            now = time.monotonic()
            self._next = max(self._next, now + delay)
            wait = self._next - now
            self._next += self.interval
        if wait > 0:
            time.sleep(wait)

    def _cache_path(self, url: str, data: bytes | None) -> Path | None:
        if self.cache is None:
            return None
        key = hashlib.sha256(url.encode() + b"\0" + (data or b"")).hexdigest()
        return self.cache / key[:2] / key

    def get(self, url: str, data: bytes | None = None, accept: str | None = None) -> bytes:
        cached = self._cache_path(url, data)
        if cached is not None and cached.exists():
            body = cached.read_bytes()
            if body == b"\0404":
                raise NotFound(url)
            return body
        headers = {"User-Agent": USER_AGENT}
        if accept:
            headers["Accept"] = accept
        if data is not None:
            headers["Content-Type"] = (
                "application/json" if data[:1] in (b"{", b"[") else "application/x-www-form-urlencoded"
            )
        delay = 0.0
        server_errors = 0
        for attempt in range(self.attempts):
            self._pace(delay)
            request = urllib.request.Request(url, data=data, headers=headers)
            try:
                with urllib.request.urlopen(request, timeout=180) as response:
                    body = response.read()
                self._store(cached, body)
                return body
            except urllib.error.HTTPError as error:
                if error.code == 404:
                    self._store(cached, b"\0404")
                    raise NotFound(url) from None
                if error.code not in RETRY_STATUS or attempt + 1 == self.attempts:
                    raise
                if error.code == 429:
                    # Rate limited: honour Retry-After, otherwise back off long.
                    retry_after = error.headers.get("Retry-After", "")
                    delay = float(retry_after) if retry_after.isdigit() else min(300.0, 15.0 * 2**attempt)
                else:
                    # A server error that persists is the resource's, not the load's.
                    server_errors += 1
                    if server_errors > 2:
                        raise
                    delay = 5.0 * server_errors
            except (urllib.error.URLError, TimeoutError, ConnectionError):
                if attempt + 1 == self.attempts:
                    raise
                delay = min(120.0, 5.0 * 2**attempt)
        raise AssertionError("unreachable")

    def _store(self, path: Path | None, body: bytes) -> None:
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            partial = path.with_suffix(".part")
            partial.write_bytes(body)
            partial.replace(path)

    def json(self, url: str, data: dict | bytes | None = None) -> dict:
        body = json.dumps(data).encode() if isinstance(data, dict) else data
        return json.loads(self.get(url, body, accept="application/json"))


def rng(seed: int, stratum: str) -> random.Random:
    """Independent, reproducible random stream per stratum."""
    return random.Random(f"{seed}:{stratum}")
