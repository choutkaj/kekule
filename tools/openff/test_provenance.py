"""Offline checks that shipped OpenFF data matches its recorded provenance."""
import hashlib
import json
import unittest
import zlib

from common import HERE, ROOT, digest, verify_sources
from package_ash import unshuffle


class ShippedDataProvenance(unittest.TestCase):
    def test_upstream_sources_match_their_lock(self):
        lock = verify_sources()
        self.assertEqual(len(lock["sources"]), 5)

    def test_normalizations_match_the_recorded_export(self):
        self.assertEqual(
            digest(ROOT / "crates/kekule-openff/data/normalizations.json"),
            "06dcbb613e5cfe75441c09b8b3b6fdc10f47b78430928a0ac835487d3d3be22c",
        )

    def test_ash_bundle_decodes_to_the_locked_model(self):
        lock = json.loads((HERE / "models.lock.json").read_text(encoding="utf-8"))
        ash = next(m for m in lock["models"] if m["name"] == "openff-gnn-am1bcc-1.0.0.pt")
        data = ROOT / "crates/kekule-openff-ash/data"
        self.assertEqual(digest(data / "model.json"), ash["bundle_sha256"]["model.json"])
        weights = unshuffle(zlib.decompress((data / "weights.planes.zlib").read_bytes()))
        self.assertEqual(hashlib.sha256(weights).hexdigest(), ash["bundle_sha256"]["weights.bin"])

    def test_locked_artifacts_match(self):
        lock = json.loads((HERE / "models.lock.json").read_text(encoding="utf-8"))
        for artifact in lock["artifacts"]:
            self.assertEqual(digest(ROOT / artifact["path"]), artifact["sha256"], artifact["path"])


if __name__ == "__main__":
    unittest.main()
