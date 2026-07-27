# File-backed trajectory codec implementation plan

> **Status:** Active implementation plan for production fixed-topology trajectory
> I/O. `ARCHITECTURE.md` remains authoritative. After the milestone is fully
> implemented and canonical feature contracts describe the shipped behavior,
> remove this plan rather than retaining it as historical documentation.

## Purpose

Implement production-quality file-backed trajectory readers and writers owned by
Molecular. The immediate consumer is MolStudio, but the capability must be a
general molecular-data I/O layer rather than a viewer-specific adapter.

The existing `molecular::trajectory` module already owns the correct semantic
model:

```text
Topology
  + ordered TrajectoryFrame state
  + reusable FrameBuffer
  + sequential reader / seekable reader / writer traits
```

This milestone supplies real file codecs behind those contracts without moving
playback, caching, rendering, or project-persistence policy into Molecular.

The design takes architectural inspiration from:

- Chemfiles: one format-agnostic trajectory interface over format-specific
  implementations and explicit format capabilities;
- MDTraj: whole-file, single-frame, and chunked/streaming access as separate
  workflows, with externally supplied topology for topology-free formats;
- GROMACS tools: explicit conversion, precision, time, periodic-box, velocity,
  force, and atom-subset semantics;
- existing Molecular invariants: exact topology identity, complete dense arrays,
  explicit units, reusable buffers, structured errors, and no silent loss.

These projects are behavioral and interoperability references. Do not copy
licensed implementation code. Test fixtures must have explicit provenance and
compatible redistribution terms.

## Architectural decision

Keep the foundational `molecular` crate dependency-light. Production codecs
belong in a new workspace companion crate:

```text
crates/molecular-trajectory-io
package: molecular-trajectory-io
Rust path: molecular_trajectory_io
```

Dependency direction is one-way:

```text
molecular <- molecular-trajectory-io <- applications such as MolStudio
```

`molecular-trajectory-io` consumes and implements the public frame, buffer,
reader, seekable-reader, writer, topology, geometry, and unit contracts from
`molecular`. The foundational crate must never depend back on the codec crate.

The codec crate owns:

- file/stream opening and format detection;
- format-specific header, frame, index, and writer implementations;
- native-unit and native-precision conversion;
- file metadata and open reports;
- checked frame-offset indices;
- path-based atomic writer convenience;
- format-specific compatibility options and diagnostics.

The foundational `molecular::trajectory` module continues to own:

- `TrajectoryFrame`, `TrajectoryFrameView`, and `FrameBuffer`;
- positions, cell, velocities, forces, time, step, observations, and properties;
- exact topology binding and atom-order assertions;
- streaming traits and generic trajectory errors;
- in-memory trajectories and memory adapters;
- buffer replacement/remapping semantics.

Do not add a second frame type, topology type, coordinate array, unit system, or
trajectory domain model in the codec crate.

## Runtime dependency policy

The production path should be pure Rust and must not require users to install
Chemfiles, GROMACS, NetCDF, a C/C++ compiler, CMake, or another system library.

Third-party pure-Rust codec dependencies are allowed only after a recorded audit
of:

- license compatibility;
- Rust-version compatibility;
- panic behavior on untrusted input;
- `unsafe` usage and soundness surface;
- resource limits and checked arithmetic;
- API stability and maintenance state;
- fixture coverage against independent implementations.

For XTC, prefer a thin private adapter over an audited compatible release of the
pure-Rust `molly` crate rather than reimplementing compressed-coordinate coding
without need. Its public types must not leak into Molecular APIs. Wrap only APIs
that can satisfy Molecular's no-panic, bounded-resource, transactional-buffer,
precision, and error contracts. If the audit finds a blocker, contribute the
required correction upstream or isolate a clean implementation behind the same
codec-private boundary.

Chemfiles and MDTraj remain reference/benchmark tools, not production runtime
dependencies. A future optional Chemfiles adapter may exist in a separately named
adapter crate, but it is not the canonical implementation and must not be needed
for ordinary Molecular or MolStudio builds.

Every new Molecular-owned crate must forbid `unsafe` code in its own source.

## Feature contracts

Create one umbrella feature and one feature per implemented format:

```text
io.trajectory.file
io.trajectory.xyz
io.trajectory.dcd
io.trajectory.xtc
io.trajectory.trr
```

Update these existing contracts where their behavior or public dependencies
change:

```text
model.trajectory
api.public-facade  # only if its actual public contract changes
```

Each codec is tracked independently. A format may be marked `experimental` while
read support, write support, variant coverage, or interoperability evidence is
incomplete. It becomes `supported` only when its documented compatibility matrix,
read/write behavior, malformed-input handling, and required tests all pass.

Do not claim support for a format merely because an enum variant or extension
mapping exists.

## Initial production scope

The first complete production milestone covers:

1. multi-frame XYZ;
2. common CHARMM/NAMD/OpenMM-style DCD;
3. GROMACS XTC;
4. GROMACS TRR.

This set gives Molecular:

- a simple text codec proving the general framework;
- a widely used fixed-record Cartesian trajectory;
- a compact lossy production trajectory;
- a full-precision portable trajectory carrying optional velocities and forces.

Explicit follow-on formats, not part of the first supported claim:

- Amber NetCDF;
- multi-model PDB;
- GRO/G96 trajectory streams;
- Amber ASCII trajectory;
- LAMMPS dump;
- TNG/H5MD and other container formats;
- transparent gzip/xz/bzip2 wrappers.

These formats must return a structured unsupported-format result until their own
feature contracts are implemented. Do not route them silently through Chemfiles
or label them supported based only on extension recognition.

## Public API direction

Exact names may be refined during implementation, but the semantics should
remain equivalent to this shape:

```rust
use molecular::trajectory::{
    AtomOrderAssertion, FrameBuffer, SeekableTrajectoryReader,
    TrajectoryFormat, TrajectoryReader, TrajectoryWriter,
};
use molecular::topology::Topology;
use molecular_trajectory_io::{
    create_trajectory_writer, open_indexed_trajectory,
    open_trajectory, FileTrajectoryMetadata, IndexedFileTrajectoryReader,
    SequentialFileTrajectoryReader, TrajectoryOpenOptions,
    TrajectoryOpenReport, TrajectoryTopologyBinding,
    TrajectoryWriteOptions,
};
```

Conceptually:

```rust
pub fn open_trajectory(
    path: impl AsRef<Path>,
    binding: TrajectoryTopologyBinding,
    options: TrajectoryOpenOptions,
) -> Result<(SequentialFileTrajectoryReader, TrajectoryOpenReport), TrajectoryError>;

pub fn open_indexed_trajectory(
    path: impl AsRef<Path>,
    binding: TrajectoryTopologyBinding,
    options: TrajectoryOpenOptions,
) -> Result<(IndexedFileTrajectoryReader, TrajectoryOpenReport), TrajectoryError>;

pub fn create_trajectory_writer(
    path: impl AsRef<Path>,
    topology: Topology,
    options: TrajectoryWriteOptions,
) -> Result<FileTrajectoryWriter, TrajectoryError>;
```

The sequential reader implements `TrajectoryReader`. The indexed reader
implements both `TrajectoryReader` and `SeekableTrajectoryReader`. The writer
implements `TrajectoryWriter` and exposes an explicit consuming `finish()`.

Format-specific modules may additionally expose concrete readers, writers, and
options:

```text
molecular_trajectory_io::xyz
molecular_trajectory_io::dcd
molecular_trajectory_io::xtc
molecular_trajectory_io::trr
```

The format-agnostic factory is a convenience over those implementations, not a
stringly typed plugin registry.

Do not add trajectory codecs to the broad Molecular prelude.

## Format identity and detection

Introduce a focused non-exhaustive `TrajectoryFormat` vocabulary, preferably in
`molecular::trajectory` because errors, reports, project descriptors, and codec
implementations share it:

```rust
#[non_exhaustive]
pub enum TrajectoryFormat {
    Xyz,
    Dcd,
    Xtc,
    Trr,
}
```

Opening supports:

```rust
pub enum TrajectoryFormatHint {
    Auto,
    Explicit(TrajectoryFormat),
}
```

`Auto` detection must not trust the filename extension alone. Use bounded
inspection of format signatures/header structure, with the extension as one
piece of evidence. Required behavior:

- explicit format bypasses extension dispatch but still validates the header;
- matching extension and signature succeed;
- known extension with conflicting signature returns a structured mismatch by
  default;
- missing extension may succeed when the signature is conclusive;
- ambiguous or insufficient evidence returns a structured detection error;
- detection reads only a bounded prefix and restores the stream position;
- unsupported compressed wrappers are reported, not misparsed as the inner
  format.

The open report records the selected format and the evidence used.

## Path and stream APIs

Path-backed files are the primary milestone, but codec implementations should
separate byte-stream logic from path convenience:

- sequential format readers should work over `Read` where practical;
- indexed readers require `Read + Seek`;
- writers require `Write` or `Write + Seek` according to format finalization;
- path factories own `BufReader<File>`/`BufWriter<File>` and attach path context
  to errors;
- explicit-format stream constructors support tests and future embedded archive
  sources without pretending an arbitrary stream has an extension.

Do not require an application to copy an embedded trajectory to a second full
in-memory buffer merely to decode it. Transparent archive integration may follow,
but the codec boundary must not make it impossible.

Async I/O is not part of this layer. Applications may run blocking readers on
worker threads. Codec types need not promise `Sync`; document `Send` behavior
only when actually established.

## Topology and atom-order binding

DCD, XTC, and TRR do not establish Molecular chemical identity. XYZ normally
provides element labels but still does not identify repeated atoms uniquely.
Therefore every initial codec requires an externally supplied exact Molecular
`Topology` and an explicit atom-order contract.

Define a checked binding conceptually equivalent to:

```rust
pub struct TrajectoryTopologyBinding {
    topology: Topology,
    atom_order: AtomOrderAssertion,
}
```

Required behavior:

- the assertion must belong to the same exact topology identity;
- the file atom count must match the topology atom count;
- the file's frame arrays are interpreted in authoritative
  `TopologyAtomIndex` order only after the assertion is validated;
- atom count alone never proves atom correspondence;
- all stronger metadata present in the format is validated;
- a metadata mismatch is an error, not a warning that silently continues.

Extend `AtomOrderAssertion` only as needed with explicit, well-named constructors
and identity queries. Suitable workflows include:

```rust
AtomOrderAssertion::from_semantic_order(&topology, &[InstanceAtomId])
AtomOrderAssertion::assert_file_uses_topology_order(&topology)
```

The second operation is deliberately explicit: it records a caller assertion for
formats that cannot prove order internally. Avoid names such as `unchecked` or a
hidden default that makes equal atom count sufficient.

For XYZ, validate the element sequence against topology in addition to the
explicit order assertion. If atom names or other reliable labels are supported by
a future format, validate them and record that evidence in the open report.

The open report should distinguish:

```text
explicit semantic-order proof
explicit dense-order assertion
element-sequence validation
format metadata unavailable
```

## File metadata and capabilities

Opening returns immutable metadata/report values rather than forcing callers to
infer capabilities from the format enum.

Conceptually:

```rust
pub struct FileTrajectoryMetadata {
    pub format: TrajectoryFormat,
    pub atom_count: usize,
    pub declared_frame_count: Option<u64>,
    pub indexed_frame_count: Option<u64>,
    pub fields: TrajectoryFieldAvailability,
    pub coordinate_encoding: CoordinateEncoding,
    pub random_access: RandomAccessCapability,
}
```

Metadata should expose, where known:

- positions availability;
- periodic-cell availability;
- velocity and force availability;
- time and step availability;
- native scalar precision;
- lossy versus lossless coordinate encoding;
- declared and verified/indexed frame counts separately;
- whether frame field presence is fixed or can vary;
- seek/index capability;
- format variant and endianness where relevant.

A header-declared frame count is not automatically a verified frame count.
Interrupted or concatenated files can disagree with headers. Indexed opening
validates the complete frame sequence and publishes a verified count. Sequential
opening may expose only a declaration until clean EOF.

`TrajectoryOpenReport` contains non-fatal facts such as an extension/signature
agreement, interpreted DCD time convention, ignored optional legacy blocks, or
lossy XTC precision. Unsupported semantics that would change scientific meaning
remain errors rather than warnings.

## Sequential and indexed access

Sequential opening must be fast and must not scan or materialize the complete
trajectory. It:

- validates the header and first-frame structural contract as needed;
- retains one open handle;
- decodes successive frames into a caller-owned `FrameBuffer`;
- maintains codec scratch storage for reuse;
- reports clean EOF only between complete frames.

Indexed opening performs a bounded full structural scan and builds an immutable
frame index containing checked byte offsets and minimal frame metadata. It must
not materialize coordinate arrays during indexing when the format allows records
to be skipped safely.

Required index behavior:

- checked `u64` file offsets;
- checked frame-count and index-size limits;
- deterministic frame order;
- validation of atom count and record structure for every frame;
- clean distinction between a valid final frame and truncation;
- exact indexed frame count;
- O(number of frames) index construction;
- O(1) offset lookup plus one-frame decode for ordinary random access;
- no reopening of the file for every frame or chunk.

`read_frame(index, destination)` must document cursor behavior. Preferred
contract: random reads do not change the sequential cursor. If a codec cannot
provide that cleanly, use a separate indexed handle rather than surprising
callers.

Persisted sidecar indices are a future optimization. Any persisted index must be
bound to a strong content fingerprint and validated before use; file path and
size alone are insufficient.

## Frame-buffer publication and transactionality

A decoder must never publish a partially decoded frame.

Add a focused core operation, conceptually:

```rust
pub struct FrameBufferData<'a> {
    pub positions: Quantity<&'a [Point3]>,
    pub cell: Option<PeriodicCell>,
    pub velocities: Option<Quantity<&'a [Vector3]>>,
    pub forces: Option<Quantity<&'a [Vector3]>>,
    pub time: Option<Quantity<f64>>,
    pub step: Option<u64>,
}

impl FrameBuffer {
    pub fn replace_from_data(
        &mut self,
        data: FrameBufferData<'_>,
    ) -> Result<(), FrameError>;
}
```

Exact names may differ, but the operation must:

1. validate topology identity, lengths, units, finite values, cell, time, and
   optional arrays before destination-visible mutation;
2. preserve the destination unchanged on any failure;
3. copy/convert into existing position, velocity, and force allocations;
4. clear stale optional fields when absent in the new frame;
5. avoid constructing an owned `TrajectoryFrame` or `Model`;
6. avoid per-frame heap allocation after reader and buffer warm-up.

File codecs decode into reader-owned reusable scratch arrays, validate a complete
frame, then publish through this operation.

EOF contract:

- clean EOF before the next frame returns `Ok(false)` and leaves the destination
  unchanged;
- EOF after any part of a frame begins is a truncation error and leaves the
  destination unchanged;
- `read_frame` out of range returns the existing structured range error and
  leaves the destination unchanged.

The existing `FrameBuffer::reset_dynamic_state` invariant remains: absent fields
from a later frame cannot leak values from an earlier frame.

## Error model

Retain one caller-facing `TrajectoryError` family, expanded with typed file and
codec context as needed. The foundational crate may define the generic context;
format implementations populate it without creating a dependency cycle.

Expected typed concepts include:

```text
TrajectoryIoOperation
- Detect
- Open
- Index
- ReadHeader
- ReadFrame
- WriteHeader
- WriteFrame
- Finish

TrajectoryCodecErrorKind
- UnknownFormat
- FormatMismatch
- InvalidHeader
- UnsupportedVariant
- TruncatedRecord
- InvalidRecordLength
- RecordMarkerMismatch
- InvalidFrame
- InconsistentAtomCount
- InconsistentMetadata
- InvalidPrecision
- ResourceLimitExceeded
- UnsupportedField
- NegativeOrUnrepresentableStep
```

Errors should carry where applicable:

- path or caller-supplied source label;
- trajectory format;
- operation;
- frame index;
- byte offset;
- expected and actual counts or record sizes;
- `std::io::ErrorKind` plus source message;
- nested `FrameError`/unit error;
- format-specific typed detail.

Do not collapse malformed data into `String`, panic on bad headers, or return
clean EOF for truncated data. Errors from untrusted input must not expose
unchecked allocation sizes or indexing panics.

Writer errors distinguish an unsupported field from a malformed field. Writers
must not silently omit velocities, forces, cell, time, step, or metadata that the
selected format cannot represent.

## Resource limits and untrusted input

Every open operation receives explicit or default-bounded `TrajectoryIoLimits`.
At minimum bound:

- atom count;
- declared/indexed frame count;
- frame byte size;
- index entry count and total index bytes;
- text line and comment/title length;
- record/block byte size;
- scratch allocation size;
- total bytes read during bounded format detection.

Requirements:

- checked arithmetic precedes allocation, multiplication, offset addition, and
  scalar-count conversion;
- file offsets use `u64` and convert to platform `usize` only after bounds checks;
- zero-atom frames and changing atom counts are rejected for ordinary Molecular
  fixed-topology trajectories;
- non-finite coordinates, vectors, cells, and times are rejected before
  publication;
- invalid enum/flag/precision values produce structured errors;
- no recursive parsing or attacker-controlled stack growth;
- no `unsafe` in Molecular-owned codec code;
- fuzz targets exercise detection, headers, indices, and frame parsing.

Default limits should permit realistic biomolecular trajectories while
preventing absurd header values from allocating attacker-controlled memory.
Document defaults and allow deliberate caller overrides.

## Unit and precision contract

All native values cross into Molecular through explicit units and are converted
once at the frame-buffer publication boundary.

Required principles:

- no implicit assumption that all formats use ångström;
- no renderer-oriented `f32` conversion in trajectory I/O;
- core coordinates remain Molecular's canonical `f64` model-length values;
- native precision and lossy encoding remain visible in metadata;
- writer precision is an explicit scientific option, not a hidden optimization;
- conversions are tested against independent reference files.

Format conventions for the initial milestone:

- XYZ: explicit read/write length-unit option; default ångström is documented in
  the open report because standard XYZ does not carry a reliable unit tag;
- DCD: coordinates are interpreted according to the implemented compatibility
  profile, normally ångström; time is not fabricated from ambiguous header data;
- XTC: positions and box are native nanometres, time is picoseconds, coordinates
  are lossy with explicit quantization precision;
- TRR: positions/box, velocities, forces, time, and scalar precision follow the
  tested GROMACS/XDR conventions and convert to Molecular units.

## Writer contract

Path-based writing is transactional by default:

1. create a temporary sibling file;
2. write header and frames;
3. explicitly `finish()` and finalize counts/headers;
4. flush and synchronize as practical;
5. atomically replace the destination where the platform permits;
6. preserve the previous destination on failure.

Dropping an unfinished writer must never publish a file as a successfully
completed trajectory. It may clean up temporary data best-effort.

Writer options include:

```text
format
overwrite policy
format-specific precision
format-specific variant
strict supported-field policy
optional declared frame count when required
```

The default supported-field policy is strict rejection. Conversion tools may
explicitly request a documented loss policy at a higher layer, but the codec
writer itself must never discard fields silently.

All frames must match the writer's exact topology identity and atom order.
Writers validate fixed atom count, units, finite values, and per-format field
consistency before changing a previously valid output.

Appending to existing trajectories, in-place repair, concatenation, and reactive
changing-topology output are not part of the first milestone.

## XYZ contract

Implement a strict multi-frame XYZ vertical slice proving the complete framework.

Reader requirements:

- bounded line parsing;
- checked nonzero atom count per frame;
- constant atom count across frames;
- element symbol validation against topology order;
- finite Cartesian coordinates;
- explicit/default length-unit policy;
- clean handling of final newline and ordinary comments;
- no invented time, step, cell, velocity, or force semantics;
- optional preservation of the bounded frame comment as namespaced frame
  metadata if the property contract supports it losslessly.

Writer requirements:

- positions and element symbols in topology order;
- explicit decimal precision;
- deterministic locale-independent formatting;
- bounded/generated comment line;
- rejection of unsupported dynamic fields by default;
- multi-frame round trip.

Extended XYZ schemas are a separate follow-on contract unless explicitly added
and tested. Do not infer cell or properties from arbitrary comment strings.

## DCD contract

Support a deliberately documented compatibility profile rather than claiming
all historical DCD dialects.

Required read coverage for the first supported profile:

- common CHARMM/NAMD/OpenMM coordinate trajectories;
- little- and big-endian 32-bit Fortran record markers;
- `CORD` trajectory headers;
- title and atom-count records with resource limits;
- ordinary all-atom frames;
- common fixed-atom/free-atom trajectories;
- common optional CHARMM/NAMD unit-cell records;
- header-declared versus indexed frame-count validation;
- truncated and mismatched record-marker detection.

Required write coverage:

- one canonical little-endian compatibility profile;
- positions and optional periodic cell;
- deterministic headers/titles;
- explicit start step/save interval/time convention;
- finalized frame count;
- no emission of fixed-atom optimization unless separately implemented and
  tested.

DCD time metadata is historically ambiguous across producers. The default reader
must preserve step information that is defensible and leave time absent unless a
format profile or explicit `DcdTimePolicy` establishes the conversion. The open
report records the interpretation. Do not silently label ambiguous values as
picoseconds.

Unit-cell angle/cosine variants and unsupported legacy records require explicit
variant detection and tests. Unknown variants return `UnsupportedVariant`, not
best-effort geometry.

The feature document must list the exact tested producer/version matrix.

## XTC contract

Implement XTC through an audited private codec adapter, preferably `molly`, while
retaining Molecular-owned semantics and errors.

Required behavior:

- validate supported XTC magic/version variants;
- fixed atom count across frames;
- positions, triclinic box, simulation step, and time;
- nanometre/picosecond conversion;
- exposure of lossy coordinate precision in metadata;
- reusable decode scratch;
- sequential read without building a full index;
- indexed offsets for random access;
- explicit writer coordinate resolution/precision;
- round trips within the declared quantization error;
- rejection of negative or unrepresentable steps;
- no panic from malformed compressed streams or small/large atom special cases.

Expose writer precision in physical terms, preferably as a positive finite length
resolution, while translating privately to the codec's native scale factor.
Document the default and test the resulting maximum error.

If the dependency exposes unbounded convenience methods or panic-prone methods,
do not call them on untrusted input. Use bounded low-level APIs or contribute a
safe upstream path.

## TRR contract

Implement a pure-Rust bounded XDR reader/writer for the tested GROMACS TRR
profile.

Required read behavior:

- XDR framing and magic/version validation;
- both single- and double-precision scalar payloads;
- fixed atom count;
- positions and triclinic box;
- optional velocities and forces;
- time and nonnegative step;
- per-frame field-size consistency;
- safe skipping/reporting of recognized optional blocks not represented by the
  initial Molecular frame contract;
- truncation and oversized-block rejection;
- indexed random access by scanned frame offsets.

Required write behavior:

- explicit single- or double-precision option;
- positions, optional box, velocities, forces, time, and step;
- stable field-presence policy across frames where required;
- explicit handling of GROMACS lambda or other supported scalar frame metadata;
- rejection of unsupported fields rather than silent loss.

If lambda is retained through `PropMap`, define and export a stable namespaced
property key and test read/write preservation. If it is promoted to a typed core
field, that must be a separately justified core API change. Do not drop a
non-default lambda silently.

## Deliberate non-goals

The first milestone does not include:

- changing-topology/reactive trajectories;
- automatic topology inference from coordinate-only files;
- atom correspondence from equal count alone;
- symmetry-aware atom reordering;
- atom-subset decoding into incomplete dense arrays;
- automatic periodic unwrapping, imaging, centering, or molecule reconstruction;
- concatenation, overlap removal, or time rewriting;
- trajectory analysis algorithms;
- playback, interpolation, prefetch, LRU caches, or GPU upload rings;
- project archive persistence or external-resource fingerprints;
- network/HTTP range readers;
- asynchronous runtime integration;
- transparent use of external Chemfiles installations;
- a generic dynamic plugin ABI.

Frame stride, time-range selection, atom subsets, concatenation, and conversion
are higher-level adapters over correct reader/writer primitives. Atom subsets in
particular require an explicit target topology and mapping; they must not weaken
complete `FrameBuffer` state.

## Tests and interoperability evidence

### Focused unit tests

For every codec, test:

- valid smallest practical file;
- multiple frames;
- exact topology/buffer identity checks;
- atom-count and available metadata mismatches;
- clean EOF versus truncated frame;
- malformed header and record sizes;
- checked overflow/resource-limit paths;
- non-finite values;
- absent optional fields clearing stale buffer state;
- destination unchanged after every late failure;
- stable position/vector pointers and capacities after warm-up;
- sequential and indexed reads producing identical frames;
- out-of-range random access;
- writer unsupported-field rejection;
- finish/atomic-output failure behavior;
- read/write round trip at the format's declared precision.

### Format-specific variants

Test at least:

- XYZ comments, whitespace, line endings, symbols, and configured units;
- DCD endianness, cells, fixed atoms, header count disagreement, record-marker
  corruption, and tested producer variants;
- XTC supported magic variants, small and large atom counts, triclinic cells,
  step/time, precision, and corrupted compressed blocks;
- TRR f32/f64, each optional vector field combination, triclinic cells, time,
  step, lambda policy, and oversized/truncated XDR blocks.

### Provenance-pinned fixtures

Store a deliberately small redistributable interoperability corpus with a
manifest recording:

```text
format and variant
generator and exact version
generator command/source script
atom count and frame count
fields and native units
expected values or hashes
license/provenance
```

Generate independent fixtures using appropriate released tools such as GROMACS,
OpenMM/NAMD-compatible DCD writers, MDTraj, and Chemfiles. Do not make those tools
runtime dependencies. Do not copy LGPL implementation code from MDTraj.

Cross-read fixtures from at least two independent implementations where practical.
Writer outputs should be accepted by an independent reader. Compare:

- atom/frame counts;
- positions within format precision;
- full triclinic cell matrices;
- velocities and forces with unit conversion;
- time and step;
- field-presence semantics;
- XTC quantization bounds;
- DCD variant interpretation.

External benchmark observations are informational and follow repository benchmark
rules. Tiny committed interoperability fixtures remain ordinary regression tests,
not optional broad-corpus benchmarks.

### Fuzzing

Add bounded fuzz targets for:

- format detection;
- XYZ frame parser;
- DCD header/record/index scanner;
- TRR/XDR frame scanner;
- Molecular-owned wrapper validation around XTC decoding where feasible.

Seed with valid and minimally corrupted fixtures. Fuzz targets must enforce small
resource limits and must not require external native libraries.

## Performance requirements

The milestone must benchmark and preserve these invariants:

- sequential open does not scan the full file;
- indexed open is O(file size) and stores O(frame count) compact offsets;
- steady-state sequential decode performs no per-frame heap allocation after
  warm-up, except where a documented third-party codec makes this unavoidable;
- readers retain one file handle rather than reopening for each frame/chunk;
- full coordinates are decoded directly into reusable scratch and then the
  caller's `FrameBuffer`;
- no `TrajectoryFrame`, `Model`, or complete in-memory `Trajectory` is constructed
  per read;
- random access decodes one requested frame after offset lookup;
- no I/O-layer conversion to renderer `f32`;
- writer buffering is bounded and `finish()` does not reread the full output.

Record benchmarks for representative small and large atom counts:

```text
sequential frames/s and MB/s
index construction time and memory
random-frame latency
allocations after warm-up where measurable
writer throughput
```

Compare against Chemfiles/MDTraj or native tools only as informational reference.
Do not distort APIs or correctness to win a benchmark.

## MolStudio integration contract

MolStudio consumes Molecular file readers directly. Molecular owns decoding,
format metadata, topology/order validation, units, cells, velocities, forces,
time, step, and structured codec errors.

MolStudio continues to own:

- `TrajectoryAsset` project IDs and provenance links;
- external versus embedded project storage policy;
- content fingerprints and relinking workflows;
- background task scheduling and cancellation;
- chunk/prefetch policy;
- LRU CPU cache;
- GPU upload ring and playback state;
- rendering-oriented `f32` packing and bounds.

After Molecular codecs satisfy the consumer contract:

- remove MolStudio's direct `chemfiles` dependency and local
  `ChemfilesTrajectorySource` ownership;
- replace local format sniffing, atom-order checks, and file errors with Molecular
  APIs;
- adapt Molecular `FrameBuffer` data to renderer chunks without re-decoding or
  reconstructing chemistry;
- retain structured capability errors only for genuinely unimplemented follow-on
  formats;
- pin MolStudio to the exact merged Molecular revision and remove any path-source
  lockfile before integration.

Required initial MolStudio consumer formats are XYZ, DCD, XTC, and TRR. Amber
NetCDF, PDB, and GRO remain visibly unavailable until their upstream feature
contracts land.

## Staged implementation strategy

This functionality is too important for one opaque unreviewable commit. Use a
draft PR with staged commits or stacked draft PRs at the boundaries below. Every
stage leaves the workspace buildable and its claimed features honestly labelled.

## Stage 0: Baseline, contracts, and dependency audit

- Read `AGENTS.md`, `ARCHITECTURE.md`, this plan, `model.trajectory`, and affected
  public contracts.
- Inventory current trajectory users in Molecular and MolStudio.
- Create feature contracts with `planned` or `experimental` status.
- Record the current four locked Molecular gates and MolStudio consumer baseline.
- Audit the proposed XTC dependency and record license/MSRV/panic/unsafe/API
  findings before adding it.
- Assemble provenance manifests for initial external fixtures.

Exit criterion: format scope, dependency direction, fixtures, and public
semantics are settled before implementation.

## Stage 1: Core publication and error contracts

- Add the minimum `molecular::trajectory` API needed by external codecs:
  transactional allocation-reusing `FrameBuffer` publication, explicit topology
  order binding helpers, format identity, and typed I/O/codec error context.
- Strengthen existing `copy_from` transactionality if the new shared kernel makes
  that practical.
- Add downstream compile tests proving a companion crate can implement the
  reader/writer traits without private access.
- Update `model.trajectory` honestly; no production codec is supported yet.

Exit criterion: file decoders can publish complete frames atomically without
owned intermediate models or duplicate frame types.

## Stage 2: Companion crate, registry, and XYZ vertical slice

- Add `molecular-trajectory-io` to the workspace.
- Implement bounded detection, path/stream constructors, metadata, reports,
  limits, sequential/indexed reader wrappers, and atomic writer scaffolding.
- Implement strict multi-frame XYZ read/write.
- Add transactionality, allocation-reuse, malformed-input, round-trip, and
  downstream public API tests.
- Mark only `io.trajectory.xyz` supported when its contract is complete.

Exit criterion: the full architecture works end to end on a transparent text
format.

## Stage 3: DCD

- Implement DCD detection, header parser, index, reader, writer, cells, fixed-atom
  compatibility, endianness, time policy, and variant diagnostics.
- Add provenance-pinned fixtures from independent producers/readers.
- Fuzz header/record scanning and benchmark sequential/indexed access.
- Mark `io.trajectory.dcd` supported only for the documented compatibility
  matrix.

Exit criterion: common production DCD files round-trip/interoperate without
silent metadata invention.

## Stage 4: TRR and shared XDR foundation

- Implement a private bounded XDR primitive layer.
- Implement TRR f32/f64 indexing, reading, writing, fields, units, and metadata.
- Add fixtures from GROMACS and an independent reader/writer.
- Fuzz XDR/frame scanning and benchmark.

Exit criterion: full-precision GROMACS trajectory state maps cleanly onto
Molecular frames.

## Stage 5: XTC

- Integrate the audited private XTC codec adapter.
- Add bounded sequential/indexed access, scratch reuse, precision metadata,
  writer resolution, error mapping, and corruption tests.
- Validate against GROMACS plus at least one independent implementation.
- Mark `io.trajectory.xtc` supported only after the documented magic/precision
  matrix passes.

Exit criterion: compact lossy trajectories are scientifically explicit and
production usable without a native library.

## Stage 6: Unified validation and MolStudio migration

- Run all locked Molecular checks, docs, package, dashboard, skills, fuzz-build,
  and focused codec tests.
- Run deliberate interoperability fixtures and record results without broad
  unsupported parity claims.
- Validate MolStudio through the untracked sibling patch.
- Replace MolStudio's local Chemfiles source with Molecular readers, remove the
  direct Chemfiles dependency, and run its locked workspace, visual, trajectory,
  and performance checks.
- Regenerate MolStudio's lockfile from the exact public Molecular revision only at
  integration.
- Remove `traj_plan.md` once canonical architecture, feature contracts, rustdoc,
  and codec docs own the final specification.

Exit criterion: MolStudio can open, index, seek, stream, and play XYZ/DCD/XTC/TRR
from a reproducibly pinned public upstream revision with no local codec substitute.

## Follow-on order

After the first milestone, recommended order is:

1. Amber NetCDF in a focused dependency/format feature;
2. GRO/G96 and multi-model PDB coordinated with Molecular structural I/O;
3. Amber ASCII and LAMMPS dump;
4. compressed stream wrappers;
5. fingerprint-bound persisted indices;
6. explicit atom-subset/topology-mapping readers;
7. concatenated multi-file readers and conversion utilities;
8. embedded project-archive trajectory blocks.

Each follow-on format receives its own contract, fixture matrix, limits, and
supported-field policy.

## Completion criteria

The initial file-backed trajectory milestone is complete only when:

- dependency direction remains acyclic and the foundational crate remains
  lightweight;
- production XYZ, DCD, XTC, and TRR contracts are individually supported;
- sequential readers do not scan whole files;
- indexed readers provide verified counts and bounded random access;
- frame publication is transactional and allocation-reusing;
- exact topology/order evidence is required;
- native units, precision, lossy encoding, and optional fields are explicit;
- malformed/truncated files return structured errors without panic;
- writers reject silent field loss and finalize atomically;
- provenance-pinned interoperability fixtures pass;
- required fuzz targets build and focused fuzz smoke is recorded where available;
- performance invariants are measured and documented;
- all required Molecular checks pass;
- MolStudio passes against the upstream implementation with its local Chemfiles
  ownership removed;
- codecs are available at a public Git revision suitable for exact pinning;
- this temporary plan has been removed in favor of canonical architecture,
  feature contracts, rustdoc, and format documentation.

## Design references

- [Chemfiles](https://chemfiles.org/) and its
  [trajectory API](https://chemfiles.org/chemfiles/latest/classes/trajectory.html)
  and [format capability matrix](https://chemfiles.org/chemfiles/latest/formats.html)
- [MDTraj trajectory loading](https://mdtraj.readthedocs.io/en/stable/load_functions.html)
  and chunked `iterload` behavior
- [GROMACS file-format reference](https://manual.gromacs.org/current/reference-manual/file-formats.html)
- [molly XTC crate documentation](https://docs.rs/molly/latest/molly/)
