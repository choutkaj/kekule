#!/usr/bin/env python3
"""RDKit reference-data generator placeholder."""

from __future__ import annotations

import argparse
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--feature", required=True)
    parser.parse_args()
    print("RDKit reference generation is not implemented for this feature yet.", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
