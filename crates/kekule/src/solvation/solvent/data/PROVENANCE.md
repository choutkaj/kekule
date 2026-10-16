# Bundled TIP3P coordinate provenance

Source: OpenMM, `wrappers/python/openmm/app/data/tip3p.pdb` at commit
`e53bdc5eab600b4137abd80e544cb100d1acac68`.

https://raw.githubusercontent.com/openmm/openmm/e53bdc5eab600b4137abd80e544cb100d1acac68/wrappers/python/openmm/app/data/tip3p.pdb

Source SHA-256: `14fd37900d627c0e258d6086a14c6084e4bec1422e9d20fccfc83c3f814fd7ee`.

This is OpenMM's pre-equilibrated 3 nm cubic box, containing 895 waters.
`tip3p.bin` retains all source coordinates without fitting, filtering, or
regeneration: 895 records of O/H1/H2, each with x/y/z signed little-endian i32
coordinates in 0.0001 nm units (36 bytes per water). Its hydrogens may extend
outside the source cell; their offsets are retained intact.

Reproduce offline with `tools/solvation/convert-tip3p.ps1 -SourcePdb PATH`.
The script checks the source checksum and atom ordering before writing.
The data is compiled into the Rust library; the converter is not a build or
runtime dependency. OpenMM attribution and license are in `LICENSE-OPENMM`.

Converted asset SHA-256: 978d946d34fbaa09a8e649d47b79689eba17d27809092b97535872810b476da2.

## Verification

OpenMM Modeller at revision `392c3d3111687867eef2d8ddae7fb3fd4470d417` supplies the algorithm conventions used here: oxygen
exclusion radius from TIP3P sigma, salt counting with 55.4 M water, and the three
padding-based lattice shapes. Kekule uses its own periodic neighbor search,
seeded ion shuffle, transactional publication, and elemental solute radii.
Its coordinates and counts are not claimed to match OpenMM byte for byte.

Algorithm source: https://github.com/openmm/openmm/blob/392c3d3111687867eef2d8ddae7fb3fd4470d417/wrappers/python/openmm/app/modeller.py
