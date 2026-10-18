# Native serialization format

Kekule's native format is the exact persistence of one canonical owner:
a `Molecule`, `Topology`, `Model`, `Ensemble`, or in-memory `Trajectory`.
Scientific formats (SMILES, Molfile, mmCIF, PDB, trajectory codecs) remain
interchange formats and are not native persistence; see
[Identity and reconstruction](../ARCHITECTURE.md#identity-and-reconstruction).
The API lives in [`kekule::native`](../crates/kekule/src/lib.rs); the codec is
in [`io/native`](../crates/kekule/src/io/native/mod.rs).

## Design

- **Exact.** A round trip preserves everything that defines molecular and
  topology layout equality: represented chemistry, definitions and instances,
  semantic IDs, dense atom order, hierarchy (label and author identifiers,
  insertion codes, component IDs, atom-site metadata, and child order), and
  classifications together with whether each was assigned explicitly. It also
  preserves properties at every scope with their units, molecular adjacency
  order, positions, cells, occupancies, B factors, ensemble weights and
  energies, frame time, step, velocities, and forces, and collection owner
  properties. Floating-point values are stored as their IEEE-754 bits, so
  coordinates come back bitwise identical (`-0.0` stays `-0.0`).
- **Identity.** Semantic IDs round-trip exactly. A decoded topology is a new
  publication and therefore gets a fresh layout identity; selections bound to
  the source topology do not bind to the decoded one.
- **Perception is optional.** The writer stores installed perception with its
  provenance (models, policies, ring-basis fallback) by default, or omits it.
  The reader installs stored perception through `Molecule::install_perception`,
  which validates it against the decoded graph, or discards it.
- **Validated reconstruction.** Decoding never trusts the input. Molecules are
  rebuilt through `MoleculeEditor` publication (connectedness, stereo, and
  normalization checks); topologies through `TopologyBuilder` publication after
  the hierarchy passes the topology's storage validation; structures through
  their checked constructors and setters. Any failure is a
  `diagnostics::FormatError` with `Format::Native`, never a panic.
- **Bounded.** Declared counts are checked against `NativeReadOptions` limits
  (reported as `LimitExceeded`) before anything is allocated, and storage grows
  only as input bytes arrive, so a short input cannot allocate much memory.
- **Streaming.** The writer and reader work on `std::io::Write` and
  `std::io::Read`. The writer encodes through a fixed 64 KiB staging buffer
  and never holds the whole file; the reader reads exactly the bytes of one
  object, so several objects can follow each other in one stream. Neither
  buffers the caller's I/O: wrap an unbuffered file in `BufReader` or
  `BufWriter`.
- **Canonical.** Equal values encode to equal bytes: owner properties and
  property columns are written in key order, and the reader rejects keys that
  are not strictly ascending.

## Versioning

The header carries one format version, currently **1**. A reader accepts only
the versions it implements and rejects any other with
`FormatErrorDetail::NativeVersion { found, supported }` of kind
`DiagnosticKind::Unsupported`. Before kekule 1.0 a version may change
incompatibly, but every incompatible change increments the version. Unknown
header flags are rejected the same way, so a newer optional section is never
silently ignored.

## Layout

All integers are little-endian and fixed width. `f64` is the little-endian
IEEE-754 bit pattern. There is no padding or alignment.

| Notation | Encoding |
| --- | --- |
| `u8`, `u16`, `u32`, `u64`, `i8`, `i32`, `i64` | Fixed-width little-endian integer |
| `f64` | 8 bytes, `f64::to_bits` little-endian |
| `bool` | `u8`, `0` or `1`; any other value is rejected |
| `count` | `u32` element count |
| `str` | `count` byte length, then that many UTF-8 bytes |
| `opt str` | `u32` length, `0xFFFF_FFFF` meaning absent, then UTF-8 bytes |
| `opt T` | `u8` tag `0` (absent) or `1` (present), then `T` |
| `bits(n)` | `ceil(n / 8)` bytes; bit `i` is bit `i % 8` of byte `i / 8`; unused high bits are zero |
| `[T; n]` | `n` consecutive `T` without a count |

### File

```text
file      := header payload trailer
header    := magic[8] version:u32 kind:u8 flags:u8 reserved:u16
magic     := 0x89 'K' 'E' 'K' 0x0D 0x0A 0x1A 0x0A
kind      := 1 molecule | 2 topology | 3 model | 4 ensemble | 5 trajectory
flags     := bit 0: perception stored; all other bits zero
reserved  := 0
trailer   := end:u32 = 0x444E_454B ("KEND") crc:u32
```

`crc` is the CRC-32 (IEEE 802.3, reflected polynomial `0xEDB88320`, as used
by zlib and PNG) of every preceding byte, header included. The decoded object
is returned only after the checksum matches
(`FormatErrorDetail::NativeChecksum` otherwise). The magic's non-ASCII first
byte and its CR LF, EOF, and LF bytes detect text-mode transfer damage.

### Molecule

Published molecule IDs are dense, so every list is in ID order and each ID is
its position.

```text
molecule   := atom_count:count atom[atom_count]
              bond_count:count bond[bond_count]
              adjacency
              element_count:count stereo_element[element_count]
              group_count:count stereo_group[group_count]
              molecule_properties
              perception?          (present when header flag 0 is set)
atom       := element:u8 isotope:u16 charge:i8 radical_electrons:u8
              radical_spin:u8 hydrogens:u8 hydrogen_count:u8 atom_map:u32
bond       := a:u32 b:u32 order:u8
adjacency  := 0                                 ascending bond-ID order
            | 1 [bond:u32; degree(atom)] for each atom in ID order
```

- `element` is the atomic number (1–118). `isotope` and `atom_map` use `0` for
  none. `radical_electrons` `0` means no radical; `radical_spin` `0` means no
  explicit spin multiplicity. `hydrogens` is `0` (inferred, count must be `0`)
  or `1` (fixed `hydrogen_count`).
- Bond `order`: `0` zero, `1` single, `2` double, `3` triple, `4` quadruple,
  `5` dative.
- `adjacency` restores each atom's incident-bond traversal order, which
  publication does not normalize. Mode `1` must list a permutation of each
  atom's incident bonds.

```text
stereo_element := 0 center:u32 carrier[4] tetrahedral:u8
                | 1 bond:u32 left:u32 right:u32 carrier carrier double_bond:u8
                | 2 axis:u32 carrier[2] axis:u8
carrier        := 0 atom:u32 | 1 (implicit hydrogen) | 2 (implicit lone pair)
stereo_group   := kind:u8 member_count:count member:u32[member_count]
```

Orientations are `0` unspecified, then `1`/`2` for clockwise/counterclockwise
(tetrahedral, axis) or together/opposite (double bond). Group kinds: `0`
absolute, `1` relative, `2` racemic, `3` and, `4` or. An element's group is
implied by group membership.

### Properties and units

```text
owner      := entry_count:count (key:str value)[entry_count]
value      := 0 bool | 1 i64 | 2 f64 unit | 3 str
table      := row_count:count column_count:count column[column_count]
column     := key:str type:u8 unit? present:bits(row_count) cell[populated]
unit       := exponent:i32[7] scale:f64 symbol:opt str
```

Keys must be strictly ascending in byte order within each owner or table.
Column `type` is `0` bool (cells `bool`), `1` int (`i64`), `2` real (`unit`
follows the type; cells `f64`), or `3` string (`str`). A stored column has
at least one present cell. Unit exponents are in the order length, mass,
time, temperature, amount, charge, angle; `scale` converts to SI. A symbol is
restored when it names one of Kekule's unit constants with the same dimension
and scale; any other symbol is dropped, which leaves the unit equal.

```text
molecule_properties  := owner table(atoms) table(bonds)
topology_properties  := owner table(instances) table(atoms) table(bonds)
                        table(chains) table(residues) table(atom_sites)
structure_properties := owner table(atoms) table(bonds)
```

### Perception

```text
perception := sections:u8 valence? rings? aromaticity? conjugation? resonance? cip?
valence    := model:u8 validation:u8 count (atom:u32 hydrogens:u8)[count]
rings      := atom_flags:bits(atom_count) bond_flags:bits(bond_count)
              basis:opt (model:u8 depth_first_fallback:bool ring_count:count ring[ring_count])
ring       := atom_count:count atom:u32[atom_count] bond_count:count bond:u32[bond_count]
aromaticity:= model:u8 count atom:u32[count] count bond:u32[count]
conjugation:= model:u8 count atom:u32[count] count bond:u32[count]
resonance  := group_count:count (count atom:u32[..] count bond:u32[..])[group_count]
cip        := count (element:u32 descriptor:u8)[count]
```

`sections` bits 0–5 mark valence, rings, aromaticity, conjugation, resonance,
and CIP in that order; other bits are rejected. Models: valence `0`
RDKit-like; validation `0` strict, `1` permissive; ring basis `0`
Figueras SSSR-like; aromaticity `0` RDKit-like, `1` MDL; conjugation `0`
RDKit-like. CIP descriptors: `0` R, `1` S, `2` r, `3` s, `4` seqTrans, `5`
seqCis, `6` E, `7` Z, `8` M, `9` P, `10` m, `11` p.

### Topology

The hierarchy is stored column by column: every value of one field for all
nodes of a kind, then the next field. Columns of strings and of ID lists keep
small reads few and fast for large systems.

```text
topology   := definition_count:count (molecule class:u8 explicit:bool)[definition_count]
              instance_count:count definition:u32[instance_count]
              atom_order
              chains residues sites
              topology_properties
atom_order := 0                                  instance order, then local atom ID
            | 1 (instance:u32 atom:u32)[atom_count]
chains     := n:count
              label:strings(n)                   required
              author:strings(n)
              residues:lists(n)                  each chain's residues in chain order
residues   := n:count
              name:strings(n)                    required
              label_comp:strings(n)
              author_comp:strings(n)
              has_label_seq:bits(n) label_seq:i32[present]
              author_seq:strings(n)
              insertion_code:strings(n)
              class:u8[n]
              explicit:bits(n)
              sites:lists(n)                     each residue's sites in residue order
sites      := n:count
              (instance:u32 atom:u32)[n]
              type_symbol:strings(n) label_asym:strings(n) auth_asym:strings(n)
              label_atom:strings(n) auth_atom:strings(n)
strings(n) := length:u32[n] bytes                lengths first, 0xFFFF_FFFF for absent;
                                                 then the present strings' UTF-8
                                                 bytes, concatenated in order
lists(n)   := length:u32[n] id:u32[total]        lengths first, then every list's
                                                 IDs, concatenated in order
```

A `required` string column has no absent entries. `label_seq` holds one `i32`
for each set bit of `has_label_seq`, in residue order.

Molecule classes: `0` protein, `1` DNA, `2` RNA, `3` carbohydrate, `4` water,
`5` ion, `6` small molecule, `7` other. Residue classes: `0` amino acid,
`1` DNA nucleotide, `2` RNA nucleotide, `3` carbohydrate, `4` water, `5` ion,
`6` other. Each residue must be listed by exactly one chain and each site by
exactly one residue; a residue's chain and a site's residue follow from those
lists. Stored classes are restored as published, explicit or inferred, without
reclassification.

### Structures

```text
model       := topology structure
ensemble    := topology owner member_count:count member[member_count]
trajectory  := topology owner frame_count:count frame[frame_count]
structure   := positions cell:opt cell occupancies b_factors structure_properties
member      := structure weight:f64 energy:opt f64
frame       := positions cell:opt cell structure_properties
               velocities:opt vectors forces:opt vectors time:opt f64 step:opt u64
positions   := (x:f64 y:f64 z:f64)[atom_count]
vectors     := (x:f64 y:f64 z:f64)[atom_count]
cell        := (x:f64 y:f64 z:f64)[3] periodic_axes:u8
occupancies := scalars
b_factors   := scalars
scalars     := 0 (absent) | 1 present:bits(atom_count) f64[populated]
```

Values are in Kekule's canonical units: positions and cell vectors in nm,
B factors in nm², velocities in nm/ps, forces in kJ/(mol·nm), energies in
kJ/mol, and time in ps. Cell vectors are the columns a, b, c;
`periodic_axes` bits 0–2 mark a, b, c. Positions and vectors have exactly
the topology's atom count; structure property tables have one row per
topology atom and bond.

## Limits

`NativeReadOptions` bounds the input bytes; the atoms and bonds of one
molecule or topology (including the instance expansion of a topology), which
also bound atom sites and residue site lists; the members or frames of one
collection, hierarchy chains and residues, property columns and owner
values; and the byte length of one string. The defaults are exported as
`DEFAULT_NATIVE_MAX_*` constants. Exceeding a limit is a `FormatError` of
kind `DiagnosticKind::Limit` carrying `LimitExceeded`. A topology can expand
few input bytes into many atoms by instantiating a large definition many
times, so set `max_atoms` to what the caller can hold when decoding untrusted
input.

## Identifier encoding

Every fixed identifier (`AtomId`, `BondId`, `StereoElementId`,
`StereoGroupId`, `MoleculeDefinitionId`, `MoleculeInstanceId`, `ChainId`,
`ResidueId`, `AtomSiteId`, `TopologyAtomIndex`, `TopologyBondIndex`) has a
stable encoding: its `raw()` `u32`, restored with `new(raw)`.
`InstanceAtomId` and `InstanceBondId` encode as `to_bits()`, a `u64` holding
the instance in the high and the local ID in the low 32 bits, restored with
`from_bits`. These encodings do not change between versions; a decoded
identifier is meaningful only against the owner it was taken from or that
owner's native round trip, whose semantic IDs are identical.
