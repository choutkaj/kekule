"""Export the two NAGL bundles that the ignored `kekule-openff` model tests read.

Run in the environment from `environment.yml`, then point the tests at the output:

    python tools/openff/export_test_bundles.py target/openff-models
    KEKULE_OPENFF_MODELS=target/openff-models cargo test -p kekule-openff --test models -- --ignored

Each bundle is checked against the fingerprints in `models.lock.json`.
"""
import argparse
import json
from pathlib import Path

from common import HERE, MODEL_HASH, ROOT, digest
from export_model import export


def export_bundles(root: Path) -> None:
    from openff.nagl_models import get_model

    lock = json.loads((HERE / "models.lock.json").read_text(encoding="utf-8"))
    rosemary = (ROOT / "crates/kekule-openff/data/rosemary.offxml").read_text(encoding="utf-8")
    for model in lock["models"]:
        name, sha = model["name"], model["checkpoint_sha256"]
        checkpoint = Path(get_model(name))
        if digest(checkpoint) != sha:
            raise ValueError(f"Checkpoint mismatch: {name}")
        directory = root / name.removesuffix(".pt")
        export(directory, checkpoint, sha, ROOT / "crates/kekule-openff/data/LICENSE-models")
        xml = rosemary.replace("openff-gnn-am1bcc-1.0.0.pt", name).replace(MODEL_HASH, sha)
        (directory / "force-field.offxml").write_text(xml, encoding="utf-8", newline="\n")
        for file, expected in model["bundle_sha256"].items():
            if digest(directory / file) != expected:
                raise ValueError(f"Bundle fingerprint mismatch: {name}/{file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("destination", type=Path, help="new directory for both bundles")
    export_bundles(parser.parse_args().destination)
