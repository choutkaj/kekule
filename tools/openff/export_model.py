"""Export a trusted NAGL checkpoint to data-only native inference files.

Run in the audited OpenFF reference environment. Python and PyTorch are exporter
dependencies only, never dependencies of the Rust engine. The checkpoint is
verified before its Python serialization is loaded.
"""
import argparse
import hashlib
import importlib.metadata
import json
import shutil
from pathlib import Path

from common import MODEL_HASH, ROOT, describe_configuration


def export(destination, checkpoint=None, checkpoint_sha256=None, license_path=None):
    import torch
    from openff.nagl import GNNModel
    from openff.nagl_models import get_model
    if checkpoint is None:
        checkpoint = get_model("openff-gnn-am1bcc-1.0.0.pt")
        checkpoint_sha256 = MODEL_HASH
        license_path = license_path or ROOT / "crates/kekule-openff/data/LICENSE-models"
    checkpoint = Path(checkpoint)
    if not checkpoint_sha256 or hashlib.sha256(checkpoint.read_bytes()).hexdigest() != checkpoint_sha256:
        raise ValueError("Checkpoint SHA-256 must be supplied and match before deserialization")
    if importlib.metadata.version("openff-nagl") != "0.6.1":
        raise ValueError("This exporter implements the openff-nagl-0.6.1 preparation profile")
    model = GNNModel.load(str(checkpoint))
    destination.mkdir(parents=True, exist_ok=False)
    weights = bytearray()
    tensors = {}
    for name, value in sorted(model.state_dict().items()):
        if value.dtype != torch.float32:
            raise ValueError(f"Unsupported tensor dtype for {name}: {value.dtype}; expected float32")
        array = value.detach().cpu().numpy().astype("<f4")
        tensors[name] = {"shape": list(array.shape), "offset": len(weights) // 4, "length": array.size}
        weights.extend(array.tobytes(order="C"))
    tables = {name: [dict(inchi=entry.inchi, mapped_smiles=entry.mapped_smiles,
                         charges=entry.property_value) for entry in table.properties.values()]
              for name, table in model.lookup_tables.items()}
    metadata = dict(schema=2, model=checkpoint.name, checkpoint_sha256=checkpoint_sha256,
                    preparation="openff-nagl-0.6.1",
                    weights_sha256=hashlib.sha256(weights).hexdigest(), tensors=tensors,
                    config=model.config.model_dump(), domain=model.chemical_domain.model_dump(),
                    lookup_tables=tables)
    (destination / "weights.bin").write_bytes(weights)
    (destination / "model.json").write_text(json.dumps(metadata, default=describe_configuration,
        separators=(",", ":"), allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    if license_path is not None:
        shutil.copyfile(license_path, destination / "LICENSE-model")
    (destination / "ATTRIBUTION.txt").write_text(
        f"Source checkpoint: {checkpoint.name}\nSHA-256: {checkpoint_sha256}\n"
        "Retain the source model's license and attribution when distributing this bundle.\n"
        "Converted to JSON and little-endian float32 tensors by kekule-openff.\n"
        "Weights and lookup values have not been retrained or symmetrized.\n",
        encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in metadata.items() if k not in ("lookup_tables", "config", "tensors")}, indent=2))
    print("lookup entries", {k: len(v) for k, v in tables.items()}, "weight bytes", len(weights))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--checkpoint-sha256", help="Required for a custom trusted checkpoint")
    parser.add_argument("--license", type=Path, help="Source model license to copy into the bundle")
    args = parser.parse_args()
    export(args.destination, args.checkpoint, args.checkpoint_sha256, args.license)
