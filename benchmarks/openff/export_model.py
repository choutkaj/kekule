"""Export the checksum-pinned Ash checkpoint to data-only native inference files.

Run in the audited OpenFF reference environment. Python and PyTorch are exporter
dependencies only, never dependencies of the Rust engine. The checkpoint is
verified before its Python serialization is loaded.
"""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

from audit import HERE, MODEL_HASH, describe_configuration


def export(destination, checkpoint=None):
    from openff.nagl import GNNModel
    from openff.nagl_models import get_model
    checkpoint = Path(checkpoint or get_model("openff-gnn-am1bcc-1.0.0.pt"))
    if hashlib.sha256(checkpoint.read_bytes()).hexdigest() != MODEL_HASH:
        raise ValueError("Ash checkpoint checksum mismatch")
    model = GNNModel.load(str(checkpoint))
    destination.mkdir(parents=True, exist_ok=True)
    weights = bytearray()
    tensors = {}
    for name, value in sorted(model.state_dict().items()):
        array = value.detach().cpu().numpy().astype("<f4")
        tensors[name] = {"shape": list(array.shape), "offset": len(weights) // 4, "length": array.size}
        weights.extend(array.tobytes(order="C"))
    tables = {name: [dict(inchi=entry.inchi, mapped_smiles=entry.mapped_smiles,
                         charges=entry.property_value) for entry in table.properties.values()]
              for name, table in model.lookup_tables.items()}
    metadata = dict(schema=1, model="openff-gnn-am1bcc-1.0.0.pt", checkpoint_sha256=MODEL_HASH,
                    weights_sha256=hashlib.sha256(weights).hexdigest(), tensors=tensors,
                    config=model.config.model_dump(), domain=model.chemical_domain.model_dump(),
                    lookup_tables=tables)
    (destination / "weights.bin").write_bytes(weights)
    (destination / "model.json").write_text(json.dumps(metadata, default=describe_configuration,
        separators=(",", ":"), allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    shutil.copyfile(HERE / "fixtures/LICENSE-models", destination / "LICENSE-models")
    (destination / "ATTRIBUTION.txt").write_text(
        "Ash openff-gnn-am1bcc-1.0.0, Copyright (c) 2023 Open Forcefield Group.\n"
        "Source: https://github.com/openforcefield/openff-nagl-models\n"
        "License: CC BY 4.0; see LICENSE-models.\n"
        "Converted to JSON and little-endian float32 tensors by kekule-openff.\n"
        "Weights and lookup values have not been retrained or symmetrized.\n",
        encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in metadata.items() if k not in ("lookup_tables", "config")}, indent=2))
    print("lookup entries", {k: len(v) for k, v in tables.items()}, "weight bytes", len(weights))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--checkpoint", type=Path)
    args = parser.parse_args()
    export(args.destination, args.checkpoint)
