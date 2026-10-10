# Attribution

**Ash** (`openff-gnn-am1bcc-1.0.0`), Copyright (c) 2023 Open Forcefield Group,
from [openff-nagl-models](https://github.com/openforcefield/openff-nagl-models).
Licensed under [CC BY 4.0](LICENSE-model).

- Source checkpoint: `openff-gnn-am1bcc-1.0.0.pt`,
  SHA-256 `7981e7f5b0b1e424c9e10a40d9e7606d96dcd3dd2b095cb4eeff6829f92238ee`.
- Converted with `tools/openff/export_model.py` (openff-nagl 0.6.1)
  to JSON metadata and little-endian float32 tensors, then packaged with
  `tools/openff/package_ash.py` as zlib-compressed byte planes.
- Decoded weights SHA-256:
  `a42f9e6325b22db821553287d2e96bb31a542df420fa5f9ca8887454ba4b9ea0`.
  `data/model.json` SHA-256:
  `80e13a479c6bba53cb94446b514b7cbf40edff01d28b28a0fd3b65e7cf967e89`.
- Changes: storage format only. Weights and lookup values have not been
  retrained, symmetrized, or otherwise modified.

Retain this attribution and the license when redistributing the model.
