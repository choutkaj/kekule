# Native serialization format

Kekule's native format is the exact persistence of one canonical owner:
a `Molecule`, `Topology`, `Structure`, `StructureSet`, or `Ensemble`. A `Structure` or `StructureSet`
may carry its structure interpretation reports, as a
`StructureInterpretation` or `StructureSetInterpretation`. A structure set
record may also store only the number of its geometries, kept elsewhere, as a
`StructureSetDescription`, and an ensemble record likewise with its
distribution, as an `EnsembleDescription`. Scientific formats (SMILES, Molfile, mmCIF, PDB, trajectory codecs) remain
interchange formats and are not native persistence; see
[Identity and reconstruction](../SPEC.md#identity-and-reconstruction).
The API lives in [`kekule::native`](../crates/kekule/src/lib.rs); the codec is
in [`io/native`](../crates/kekule/src/io/native/mod.rs).

## Design

- **Exact.** A round trip preserves everything that defines molecular and
  topology layout equality: represented chemistry, definitions and instances,
  semantic IDs, dense atom order, hierarchy (label and author identifiers,
  insertion codes, component IDs, atom-site metadata, and child order), and
  classifications together with whether each was assigned explicitly. It also
  preserves properties at every scope with their units, molecular adjacency
  order, positions, cells, occupancies, B factors, and a structure set's geometries
  with their velocities and forces, its order (times and steps, steps, or
  coordinates with their unit), and its owner properties, and an ensemble's
  distribution (samples, log weights, each sample's state, each state's
  representative, energy, and owner properties, and the temperature). Interpretation reports
  round-trip field by field and compare equal. Floating-point values are
  stored as their IEEE-754 bits, so coordinates come back bitwise identical
  (`-0.0` stays `-0.0`).
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

The header carries one format version, currently **6**. A reader accepts only
the versions it implements and rejects any other with
`FormatErrorDetail::NativeVersion { found, supported }` of kind
`DiagnosticKind::Unsupported`. Before kekule 1.0 a version may change
incompatibly, but every incompatible change increments the version. Unknown
header flags are rejected the same way, so a newer optional section is never
silently ignored.

The current version is 6 and reads no earlier version. Version 6 keeps the
distribution in the ensemble record, kind 11, apart from the structure set
record, and stores natural logarithms of weights.

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
kind      := 1 molecule | 2 topology | 3 structure | 7 interpreted structure
           | 9 structure set | 10 interpreted structure set | 11 ensemble
flags     := bit 0: perception stored; all other bits zero
reserved  := 0
trailer   := end:u32 = 0x444E_454B ("KEND") crc:u32
```

Kinds 4, 5, 6, and 8 are retired and, like any other unlisted kind, rejected
as unknown.

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
geometry_properties := owner table(atoms) table(bonds)
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
resonance  := model:u8 group_count:count (count atom:u32[..] count bond:u32[..])[group_count]
cip        := count (element:u32 descriptor:u8)[count]
```

`sections` bits 0–5 mark valence, rings, aromaticity, conjugation, resonance,
and CIP in that order; other bits are rejected. Models: valence `0`
RDKit-like; validation `0` strict, `1` permissive; ring basis `0`
Figueras SSSR-like; aromaticity `0` RDKit-like, `1` MDL; conjugation `0`
RDKit-like; resonance `0` charge moves. CIP descriptors: `0` R, `1` S, `2` r, `3` s, `4` seqTrans, `5`
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
structure      := topology geometry_state
structure_set  := set_head geometries
ensemble       := set_head distribution geometries
set_head       := topology owner geometry_count:count order
order          := 0                                (none)
               | 1 times:opt [f64; geometry_count] steps:opt [u64; geometry_count]
                                                   (time)
               | 2 [u64; geometry_count]           (steps)
               | 3 unit [f64; geometry_count]      (a coordinate in that unit)
distribution   := members log_weights:opt [f64; sample_count]
                  state:u32[sample_count]
                  state_count:count state[state_count]
                  temperature:opt f64
members        := 0                                (every geometry)
               | 1 sample_count:count geometry:u32[sample_count]
state          := representative:u32 energy:opt f64 owner
geometries     := 0 geometry[geometry_count]       (inline)
               | 1                                 (stored elsewhere)
geometry       := geometry_state velocities:opt vectors forces:opt vectors
geometry_state := positions cell:opt cell occupancies b_factors geometry_properties
positions      := (x:f64 y:f64 z:f64)[atom_count]
vectors        := (x:f64 y:f64 z:f64)[atom_count]
cell           := (x:f64 y:f64 z:f64)[3] periodic_axes:u8
occupancies    := scalars
b_factors      := scalars
scalars        := 0 (absent) | 1 present:bits(atom_count) f64[populated]
```

Values are in Kekule's canonical units: positions and cell vectors in nm,
B factors in nm², velocities in nm/ps, forces in kJ/(mol·nm), energies in
kJ/mol, temperatures in K, and time in ps. Cell vectors are the columns a, b, c;
`periodic_axes` bits 0–2 mark a, b, c. Positions and vectors have exactly
the topology's atom count; geometry property tables have one row per
topology atom and bond.

A structure set stores its geometries in set order, after what describes them. A
time order stores times, steps, or both, each for every geometry; a step order
a step per geometry; a coordinate order a coordinate per geometry in its unit,
stored as property units are. Times and coordinates are finite.

An ensemble is a structure set and a distribution over its geometries; a
structure set record states no weights or states. A distribution's samples
are every geometry, or the geometries that a strictly ascending list of indices
names; `sample_count` is the geometry count for every geometry. Log weights,
the natural logarithm of each sample's weight, one per sample, are finite;
without them samples weigh the same. Each sample stores the index of its
state, and each state the sample position of its representative. The reader
validates the partition: every state index names a stored state, and every
representative is a sample of its own state, so no state is empty.

The geometries section is inline, or tag `1` for geometries stored elsewhere: such
a structure set record is a `StructureSetDescription`, whose `geometry_count`
bounds its order as stored geometries would, and such an ensemble record an
`EnsembleDescription`, whose `geometry_count` bounds its distribution too.
`StructureSetDescription::of` and `EnsembleDescription::of` describe a set or
ensemble in memory or read from a trajectory file, so a host stores a file's
order and distribution beside the file. An interpreted structure set stores
its geometries inline.

### Interpretation reports

An interpreted structure is a structure followed by its report; an interpreted
structure set is a structure set followed by one report per geometry, in geometry
order. Reports are stored as they are: their IDs name atoms, chains,
residues, and molecule instances of the structure they were interpreted with,
and the reader does not resolve them.

Report columns use two further encodings. `packed(n)` holds `n` unsigned
integers in the narrowest width that holds the largest, and `dict(n)` holds
a string column through a table of its distinct values, in order of first
appearance.

```text
packed(n)  := width:u8 value[n]                  width 1, 2, 4, or 8 bytes,
                                                 little-endian values
dict(n)    := value_count:count strings(value_count) index:packed(n)
                                                 value_count <= n; every index
                                                 names a table value
list(T)    := count T[count]
usize      := u64
```

A `required` dictionary has no absent table values.

```text
report      := block_name:str selected_model:opt str
               entity_definitions:usize coordinate_models:usize
               macromolecules:usize small_molecules:usize
               solvent_molecules:usize applied_connections:usize
               connectivity_candidates:usize
               ignored_coordinate_models:list(str)
               instances
               coordination:list(coordination)
               alternate_locations:list(decision)
               issues:list(issue)
               polymer_entities:list(entity)
               assemblies:list(assembly)
instances   := n:count
               molecule:packed(n)
               coordinate_model:dict(n)          required
               asym_count:packed(n) asym_id:dict(Σ)            required
               entity_count:packed(n) entity_id:dict(Σ)        required
               kind_count:packed(n) entity_kind:kinds(Σ)
               atom_count:packed(n) bond_count:packed(n)
               atoms(Σ atom_count) bonds(Σ bond_count)
kinds(m)    := code:packed(m) other:dict(code 4 rows)       required
atoms(m)    := atom:refs(m)
               source_line:packed(m) source_order:packed(m)
               atom_site_id:strings(m)
               type_symbol:dict(m)               required
               label_atom_name:dict(m) atom_name:dict(m) auth_atom_name:dict(m)
               label_component_id:dict(m) component_id:dict(m)
               auth_component_id:dict(m)
               label_asym_id:dict(m) asym_id:dict(m) auth_asym_id:dict(m)
               entity_id:dict(m)
               residue_key:dict(m)               required
               author_sequence_id:dict(m) insertion_code:dict(m)
               selected_alternate_location:dict(m)
               entity_kind:kinds(m)
               has_label_sequence:bits(m) label_sequence:i32[present]
               has_occurrence:bits(m) occurrence:packed(present)
bonds(m)    := atoms:refs(2m) evidence:packed(m)
               template_component:dict(evidence 1 rows)      required
               template_variant:dict(evidence 1 rows)
               connection_id:dict(evidence 4 rows)
               connection_type:dict(evidence 4 rows)         required
refs(m)     := 0 local:packed(m)                 atom IDs of each row's instance
             | 1 id:packed(m)                    InstanceAtomId::to_bits
```

`atom_name`, `component_id`, and `asym_id` are required. The atoms and bonds
of every instance are stored in instance order, each instance's in its own
order. `refs` mode `0` holds local atom IDs when every reference is to the
molecule instance of its provenance; mode `1` holds full IDs.
Entity kinds: `0` polymer, `1` branched, `2` non-polymer, `3` water, `4`
another value, stored in `other`. Bond evidence: `0` component table, `1`
template, `2` polymer link, `3` branch link, `4` connection, `5` distance,
`6` hydrogen distance.

```text
coordination := a:u64 b:u64 connection_id:opt str connection_type:str
                source_line:opt usize
decision     := residue available:list(str) selected:list(str) reason
                selected_lines:count packed omitted_lines:count packed
reason       := 0 highest occupancy | 1 exact selection | 2 preferred label
              | 3 preferred label unavailable | 4 id:str (source conformation)
residue      := model_id:str asym_id:str position insertion_code:opt str
position     := 0 label:i32 | 1 author:str | 2 author:str repeat:usize
              | 3 component_id:str index:usize
entity       := id:str polymer_type:opt (code:u8 other:str when 9)
                monomers_source:u8 sequence:str sequence_source:u8
                list(monomer) list(entity_instance)
monomer      := number:i32 component_id:str alternatives:list(str)
                code:str origin:u8 parents:list(str)
entity_instance := asym_id:str author_chain_id:opt str chain:opt u32
                   list(sequence_position)
sequence_position := residue:opt u32 author_seq_id:opt str insertion_code:opt str
assembly     := id:str details:opt str oligomeric_details:opt str
                oligomeric_count:opt usize list(generator)
generator    := chains:list(str) list(operation)
operation    := operators:list(str) rotation:f64[9] translation:f64[3]
```

An operation's rotation is its three columns, each `x y z`: a proper
rotation, each entry of RᵀR and its determinant within 10⁻³ of the
identity's. Its translation is in the canonical length unit.

Coordination atoms are `InstanceAtomId::to_bits`. Polymer types: `0`
polypeptide(L), `1` polypeptide(D), `2` polydeoxyribonucleotide, `3`
polyribonucleotide, `4` DNA/RNA hybrid, `5` peptide nucleic acid, `6`
cyclic-pseudo-peptide, `7` polysaccharide(D), `8` polysaccharide(L), `9`
another value. Sequence sources: `0` stated, `1` derived. Code origins: `0`
component, `1` parent, `2` sequence, `3` unknown.

An issue is its kind, then its fields in declaration order, with lists as
`list`, line lists as `count packed`, residues as `residue`, entity kinds as
`u8` (`str` follows `4`), charges as `i8`, and partners as `u8`.

| Kind | Issue | Kind | Issue |
| --- | --- | --- | --- |
| 0 | `EntityTypeInferred` | 10 | `AlternateLocationOccupancyTie` |
| 1 | `ElementsInferredFromAtomNames` | 11 | `ConnectivityCandidatesInferred` |
| 2 | `ConnectionsIgnored` | 12 | `ResidueBondsInferred` |
| 3 | `ConnectionUnresolved` | 13 | `HydrogenAttachedByDistance` |
| 4 | `ConnectionAmbiguous` | 14 | `ChainBreak` |
| 5 | `CoordinateModelIgnored` | 15 | `ChargesFromBonds` |
| 6 | `AlternateLocationOmitted` | 16 | `ProtonationFromHydrogens` |
| 7 | `AlternateLocationIncomplete` | 17 | `ResidueFragmented` |
| 8 | `AlternateLocationOccupancyMissing` | 18 | `CellOmitted` |
| 9 | `AlternateLocationCombinationRejected` | 19 | `AssemblyOmitted` |

Connection resolution reasons: `0` unsupported symmetry (`str`), `1` missing
selector, `2` conflicting selector values (`u8`: `0` label alternate-location
aliases, `1` label and author alternate locations), `3` conflicting label and
author selectors, `4` alternate location omitted (`str`), `5` no matching
atom, `6` multiple matching atoms. A cell omission is its reason (`0`
incomplete, `1` space group with its `str` name, `2` placeholder), its line
(`usize`), and its column (`opt usize`). An assembly omission is its
assembly (`opt str`), its reason (`0` unidentified, `1` no generators, `2`
expression, `3` unknown operator, `4` unreadable operator, `5` not a
rotation, `6` too many operations; `2` to `5` with a `str`), and its line
(`opt usize`).

## Limits

`NativeReadOptions` bounds the input bytes; the atoms and bonds of one
molecule or topology (including the instance expansion of a topology), which
also bound atom sites and residue site lists; the geometries, samples, and states of one
structure set, hierarchy chains and residues, property columns and owner
values, and the entries of report lists; and the byte length of one string.
The atoms and bonds of one report's instance provenance count against the
atom and bond limits. The defaults are exported as
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
