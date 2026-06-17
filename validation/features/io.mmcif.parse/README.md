# io.mmcif.parse Validation Fixtures

These fixtures are manual PDBx/mmCIF atom-site inputs, not Biopython-generated goldens.

They focus on parser durability:

- minimal `_atom_site` loops
- label and author identifiers
- alternate locations
- insertion codes
- heterogens and waters
- multiple models
- quoted values and text blocks

The parser must keep raw parsing separate from bond perception and sanitization.
