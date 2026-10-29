//! Native decoder: one object from a `Read` source.
//!
//! The input is untrusted. Every count is checked against the read options
//! before storage is reserved, storage grows only as bytes arrive, and every
//! value reaches its owner through the owner's validating constructors, so a
//! malformed input yields a [`FormatError`], never a panic.

use std::io::{ErrorKind, Read};
use std::num::{NonZeroU16, NonZeroU32};
use std::sync::Arc;

use super::crc::Crc32;
use super::encode::{ABSENT_STRING, BASE_DIMENSIONS, END_MARKER};
use super::{NativeKind, NativeObject, NativeReadOptions, NATIVE_FORMAT_VERSION, NATIVE_MAGIC};
use crate::diagnostics::{
    DiagnosticKind, Format, FormatError, FormatErrorDetail, LimitExceeded, Resource,
};
use crate::geometry::{PeriodicCell, Point3, Vector3};
use crate::molecule::{
    Atom, AtomId, AtomRadical, AxisOrientation, AxisStereo, BondId, BondOrder,
    DoubleBondOrientation, DoubleBondStereo, Element, ImplicitHydrogens, Molecule, MoleculeEditor,
    StereoCarrier, StereoDescriptor, StereoElement, StereoElementId, StereoElementKind,
    StereoGroup, StereoGroupKind, TetrahedralOrientation, TetrahedralStereo,
};
use crate::perception::{
    aromaticity::AromaticityModel,
    conjugation::ConjugationModel,
    resonance::ResonanceGroup,
    rings::{Ring, RingBasisModel, RingBasisState, RingMembership, RingSet},
    valence::{ValenceModel, ValenceValidation},
    Perception,
};
use crate::properties::{
    MoleculeProperties, OwnerProperties, PropertyColumn, PropertyKey, PropertyValue,
    RawPropertyTable, StructureProperties, TopologyProperties,
};
use crate::structure::{
    Ensemble, EnsemblePartition, EnsembleSample, EnsembleState, Forces, Model, Positions, Series,
    SeriesAxis, SeriesItem, Trajectory, TrajectoryFrame, Velocities,
};
use crate::topology::{
    AtomSiteId, AtomSiteMetadata, Hierarchy, InstanceAtomId, MoleculeClass, MoleculeDefinitionId,
    MoleculeInstanceId, RawAtomSite, RawChain, RawDefinition, RawHierarchyError, RawResidue,
    ResidueClass, ResidueId, Topology, TopologyParts,
};
use crate::units::{
    unit_with_symbol, Dimension, Quantity, Unit, CANONICAL_ENERGY_UNIT, CANONICAL_FORCE_UNIT,
    CANONICAL_LENGTH_UNIT, CANONICAL_TEMPERATURE_UNIT, CANONICAL_TIME_UNIT,
    CANONICAL_VELOCITY_UNIT, SQUARE_NANOMETER,
};

/// Bytes read per chunk of a large array; a multiple of 24 and 8.
const CHUNK_BYTES: usize = 24 * 2730;

type Result<T> = std::result::Result<T, FormatError>;

pub(crate) fn read_native(
    reader: &mut impl Read,
    options: NativeReadOptions,
) -> Result<NativeObject> {
    let mut decoder = Decoder {
        source: reader,
        offset: 0,
        crc: Crc32::new(),
        options,
        stored_perception: false,
        scratch: Vec::new(),
    };
    let kind = decoder.header()?;
    let object = match kind {
        NativeKind::Molecule => NativeObject::Molecule(Box::new(decoder.molecule()?)),
        NativeKind::Topology => NativeObject::Topology(decoder.topology()?),
        NativeKind::Model => NativeObject::Model(decoder.model()?),
        NativeKind::Ensemble => NativeObject::Ensemble(decoder.ensemble()?),
        NativeKind::Trajectory => NativeObject::Trajectory(decoder.trajectory()?),
        NativeKind::Series => NativeObject::Series(decoder.series()?),
    };
    decoder.trailer()?;
    Ok(object)
}

struct Decoder<'r, R: Read> {
    source: &'r mut R,
    offset: u64,
    crc: Crc32,
    options: NativeReadOptions,
    stored_perception: bool,
    scratch: Vec<u8>,
}

impl<R: Read> Decoder<'_, R> {
    // Errors.

    fn error(&self, kind: DiagnosticKind, message: impl Into<String>) -> FormatError {
        FormatError::parse(Format::Native, kind, message).at_byte(self.offset)
    }

    fn rejected(&self, kind: DiagnosticKind, message: impl Into<String>) -> FormatError {
        FormatError::interpret(Format::Native, kind, message).at_byte(self.offset)
    }

    fn invalid(&self, message: impl Into<String>) -> FormatError {
        self.error(DiagnosticKind::InvalidValue, message)
    }

    fn limit(&self, resource: Resource, limit: usize, observed: usize) -> FormatError {
        self.error(
            DiagnosticKind::Limit,
            format!("native input exceeds the {resource} limit"),
        )
        .exceeding(LimitExceeded::new(resource, limit, observed))
    }

    // Primitive reads.

    fn fill(&mut self, buffer: &mut [u8]) -> Result<()> {
        let end = usize::try_from(self.offset)
            .ok()
            .and_then(|offset| offset.checked_add(buffer.len()));
        let max = self.options.max_input_bytes();
        if end.is_none_or(|end| end > max) {
            return Err(self.limit(Resource::InputBytes, max, end.unwrap_or(usize::MAX)));
        }
        self.source.read_exact(buffer).map_err(|error| {
            if error.kind() == ErrorKind::UnexpectedEof {
                self.error(DiagnosticKind::Syntax, "native input ends early")
            } else {
                self.error(DiagnosticKind::Io, "reading native input failed")
                    .with_source(error)
            }
        })?;
        self.crc.update(buffer);
        self.offset += buffer.len() as u64;
        Ok(())
    }

    fn array<const N: usize>(&mut self) -> Result<[u8; N]> {
        let mut bytes = [0; N];
        self.fill(&mut bytes)?;
        Ok(bytes)
    }

    fn u8(&mut self) -> Result<u8> {
        Ok(self.array::<1>()?[0])
    }

    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.array()?))
    }

    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.array()?))
    }

    fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.array()?))
    }

    fn i64(&mut self) -> Result<i64> {
        Ok(i64::from_le_bytes(self.array()?))
    }

    fn f64(&mut self) -> Result<f64> {
        Ok(f64::from_bits(u64::from_le_bytes(self.array()?)))
    }

    fn bool(&mut self) -> Result<bool> {
        match self.u8()? {
            0 => Ok(false),
            1 => Ok(true),
            other => Err(self.invalid(format!("invalid boolean {other}"))),
        }
    }

    fn tag(&mut self, what: &str) -> Result<bool> {
        match self.u8()? {
            0 => Ok(false),
            1 => Ok(true),
            other => Err(self.invalid(format!("invalid {what} tag {other}"))),
        }
    }

    /// A count checked against `limit` before anything is reserved for it.
    fn count(&mut self, limit: usize, resource: Resource) -> Result<usize> {
        let count = self.u32()? as usize;
        if count > limit {
            return Err(self.limit(resource, limit, count));
        }
        Ok(count)
    }

    fn string_bytes(&mut self, len: usize) -> Result<String> {
        let max = self.options.max_string_bytes();
        if len > max {
            return Err(self.limit(Resource::ValueBytes, max, len));
        }
        let mut bytes = vec![0; len];
        self.fill(&mut bytes)?;
        String::from_utf8(bytes).map_err(|_| self.invalid("a string is not valid UTF-8"))
    }

    fn str(&mut self) -> Result<String> {
        let len = self.u32()?;
        if len == ABSENT_STRING {
            return Err(self.invalid("a required string is absent"));
        }
        self.string_bytes(len as usize)
    }

    fn opt_str(&mut self) -> Result<Option<String>> {
        match self.u32()? {
            ABSENT_STRING => Ok(None),
            len => self.string_bytes(len as usize).map(Some),
        }
    }

    fn opt_f64(&mut self) -> Result<Option<f64>> {
        Ok(if self.tag("optional value")? {
            Some(self.f64()?)
        } else {
            None
        })
    }

    /// Reads `bytes` input bytes in bounded chunks, passing each to `visit`.
    fn chunks(
        &mut self,
        bytes: usize,
        mut visit: impl FnMut(&[u8]) -> std::result::Result<(), String>,
    ) -> Result<()> {
        let mut remaining = bytes;
        let mut buffer = std::mem::take(&mut self.scratch);
        while remaining > 0 {
            let len = remaining.min(CHUNK_BYTES);
            buffer.resize(len, 0);
            self.fill(&mut buffer)?;
            visit(&buffer).map_err(|message| self.invalid(message))?;
            remaining -= len;
        }
        self.scratch = buffer;
        Ok(())
    }

    fn bits(&mut self, len: usize) -> Result<Vec<bool>> {
        let mut flags = Vec::with_capacity(len);
        self.chunks(len.div_ceil(8), |bytes| {
            for &byte in bytes {
                for bit in 0..8 {
                    if flags.len() < len {
                        flags.push(byte & (1 << bit) != 0);
                    } else if byte & (1 << bit) != 0 {
                        return Err("unused bits of a bit set are not zero".into());
                    }
                }
            }
            Ok(())
        })?;
        Ok(flags)
    }

    /// A list of raw IDs whose count is at most `limit`.
    fn ids(&mut self, limit: usize, resource: Resource) -> Result<Vec<u32>> {
        let count = self.count(limit, resource)?;
        let mut ids = Vec::with_capacity(count.min(CHUNK_BYTES / 4));
        self.chunks(count.saturating_mul(4), |bytes| {
            ids.extend(
                bytes
                    .as_chunks::<4>()
                    .0
                    .iter()
                    .map(|raw| u32::from_le_bytes(*raw)),
            );
            Ok(())
        })?;
        Ok(ids)
    }

    fn u32_array(&mut self, count: usize) -> Result<Vec<u32>> {
        let mut values = Vec::with_capacity(count.min(CHUNK_BYTES));
        self.chunks(count.saturating_mul(4), |bytes| {
            values.extend(
                bytes
                    .as_chunks::<4>()
                    .0
                    .iter()
                    .map(|raw| u32::from_le_bytes(*raw)),
            );
            Ok(())
        })?;
        Ok(values)
    }

    /// A string column of `count` entries; absent entries are rejected when
    /// `required`.
    fn strings(&mut self, count: usize, required: bool) -> Result<Vec<Option<String>>> {
        let lengths = self.u32_array(count)?;
        let max = self.options.max_string_bytes();
        let mut total = 0usize;
        for &length in &lengths {
            if length == ABSENT_STRING {
                if required {
                    return Err(self.invalid("a required string is absent"));
                }
            } else if length as usize > max {
                return Err(self.limit(Resource::ValueBytes, max, length as usize));
            } else {
                total = total.saturating_add(length as usize);
            }
        }
        let mut bytes = Vec::with_capacity(total.min(CHUNK_BYTES));
        self.chunks(total, |chunk| {
            bytes.extend_from_slice(chunk);
            Ok(())
        })?;
        let mut strings = Vec::with_capacity(count);
        let mut start = 0;
        for length in lengths {
            if length == ABSENT_STRING {
                strings.push(None);
                continue;
            }
            let end = start + length as usize;
            let text = bytes
                .get(start..end)
                .and_then(|bytes| std::str::from_utf8(bytes).ok())
                .ok_or_else(|| self.invalid("a string is not valid UTF-8"))?;
            strings.push(Some(text.to_owned()));
            start = end;
        }
        Ok(strings)
    }

    /// A column of `count` ID lists whose total length is at most `limit`.
    fn lists<T>(
        &mut self,
        count: usize,
        limit: usize,
        resource: Resource,
        id: impl Fn(u32) -> T,
    ) -> Result<Vec<Vec<T>>> {
        let lengths = self.u32_array(count)?;
        let total = lengths.iter().fold(0usize, |total, &length| {
            total.saturating_add(length as usize)
        });
        if total > limit {
            return Err(self.limit(resource, limit, total));
        }
        let ids = self.u32_array(total)?;
        let mut lists = Vec::with_capacity(count);
        let mut start = 0;
        for length in lengths {
            let end = start + length as usize;
            lists.push(
                ids.get(start..end)
                    .unwrap_or_default()
                    .iter()
                    .map(|&raw| id(raw))
                    .collect(),
            );
            start = end;
        }
        Ok(lists)
    }

    // Framing.

    fn header(&mut self) -> Result<NativeKind> {
        let magic = self.array::<8>()?;
        if magic != NATIVE_MAGIC {
            return Err(self
                .error(
                    DiagnosticKind::Syntax,
                    "input is not in Kekule's native format",
                )
                .with_detail(FormatErrorDetail::NativeMagic));
        }
        let version = self.u32()?;
        if version != NATIVE_FORMAT_VERSION {
            return Err(self
                .error(
                    DiagnosticKind::Unsupported,
                    format!(
                        "native format version {version} is not supported; \
                         this reader implements version {NATIVE_FORMAT_VERSION}"
                    ),
                )
                .with_detail(FormatErrorDetail::NativeVersion {
                    found: version,
                    supported: NATIVE_FORMAT_VERSION,
                }));
        }
        let code = self.u8()?;
        let kind = NativeKind::from_code(code)
            .ok_or_else(|| self.invalid(format!("unknown native object kind {code}")))?;
        let flags = self.u8()?;
        if flags & !1 != 0 {
            return Err(self.error(
                DiagnosticKind::Unsupported,
                format!("unknown native header flags {flags:#04x}"),
            ));
        }
        self.stored_perception = flags & 1 != 0;
        let reserved = u16::from_le_bytes(self.array()?);
        if reserved != 0 {
            return Err(self.error(
                DiagnosticKind::Unsupported,
                "reserved native header bytes are not zero",
            ));
        }
        Ok(kind)
    }

    fn trailer(&mut self) -> Result<()> {
        if self.u32()? != END_MARKER {
            return Err(self.error(
                DiagnosticKind::Syntax,
                "native object is not followed by its end marker",
            ));
        }
        let computed = self.crc.finish();
        let stored = self.u32()?;
        if stored != computed {
            return Err(self
                .error(
                    DiagnosticKind::Inconsistent,
                    "native input checksum does not match its content",
                )
                .with_detail(FormatErrorDetail::NativeChecksum { stored, computed }));
        }
        Ok(())
    }

    // Molecules.

    fn molecule(&mut self) -> Result<Molecule> {
        let atom_count = self.count(self.options.max_atoms(), Resource::Atoms)?;
        let mut editor = MoleculeEditor::new();
        for _ in 0..atom_count {
            let atom = self.atom()?;
            editor.add_atom(atom).map_err(|error| {
                self.rejected(DiagnosticKind::Chemistry, "invalid atom")
                    .with_source(error)
            })?;
        }
        let bond_count = self.count(self.options.max_bonds(), Resource::Bonds)?;
        let mut degrees = vec![0usize; atom_count];
        for _ in 0..bond_count {
            let record = self.array::<9>()?;
            let a = u32::from_le_bytes([record[0], record[1], record[2], record[3]]);
            let b = u32::from_le_bytes([record[4], record[5], record[6], record[7]]);
            let order = bond_order(record[8])
                .ok_or_else(|| self.invalid(format!("unknown bond order {}", record[8])))?;
            editor
                .add_bond(AtomId::new(a), AtomId::new(b), order)
                .map_err(|error| {
                    self.rejected(DiagnosticKind::InvalidReference, "invalid bond")
                        .with_source(error)
                })?;
            // The editor accepted both endpoints, so they index `degrees`.
            for endpoint in [a, b] {
                if let Some(degree) = degrees.get_mut(endpoint as usize) {
                    *degree += 1;
                }
            }
        }
        let adjacency = if self.tag("adjacency")? {
            let mut lists = Vec::with_capacity(atom_count);
            for &degree in &degrees {
                let mut list = Vec::with_capacity(degree);
                for _ in 0..degree {
                    list.push(BondId::new(self.u32()?));
                }
                lists.push(list);
            }
            Some(lists)
        } else {
            None
        };
        let element_limit = atom_count.saturating_add(bond_count);
        let element_count = self.count(element_limit, Resource::Values)?;
        for _ in 0..element_count {
            let kind = self.stereo_element()?;
            editor
                .add_stereo_element(StereoElement::new(kind))
                .map_err(|error| {
                    self.rejected(DiagnosticKind::Chemistry, "invalid stereo element")
                        .with_source(error)
                })?;
        }
        let group_count = self.count(element_count, Resource::Values)?;
        for _ in 0..group_count {
            let code = self.u8()?;
            let kind = group_kind(code)
                .ok_or_else(|| self.invalid(format!("unknown stereo group kind {code}")))?;
            let members = self
                .ids(element_count, Resource::Values)?
                .into_iter()
                .map(StereoElementId::new)
                .collect();
            editor
                .add_stereo_group(StereoGroup { kind, members })
                .map_err(|error| {
                    self.rejected(DiagnosticKind::Chemistry, "invalid stereo group")
                        .with_source(error)
                })?;
        }
        let owner = self.owner()?;
        let atoms = self.table(atom_count)?;
        let bonds = self.table(bond_count)?;
        let perception = if self.stored_perception {
            Some(self.perception(atom_count, bond_count, element_count)?)
        } else {
            None
        };
        let mut molecule = editor.finish().map_err(|error| {
            self.rejected(
                DiagnosticKind::Chemistry,
                "stored chemistry does not publish as a molecule",
            )
            .with_source(error)
        })?;
        if let Some(lists) = adjacency {
            for (index, list) in lists.into_iter().enumerate() {
                if !molecule.restore_adjacency_order(AtomId::new(index as u32), list) {
                    return Err(self.rejected(
                        DiagnosticKind::Inconsistent,
                        format!(
                            "stored adjacency of atom {index} is not a permutation of its bonds"
                        ),
                    ));
                }
            }
        }
        let mut properties = MoleculeProperties::new(atom_count, bond_count);
        *properties.owner_mut() = owner;
        *properties.atoms_mut() = atoms;
        *properties.bonds_mut() = bonds;
        molecule.properties = properties;
        if let Some(perception) = perception.filter(|_| self.options.perception()) {
            molecule.install_perception(perception).map_err(|error| {
                self.rejected(
                    DiagnosticKind::Inconsistent,
                    "stored perception does not fit the stored chemistry",
                )
                .with_source(error)
            })?;
        }
        Ok(molecule)
    }

    fn atom(&mut self) -> Result<Atom> {
        let record = self.array::<12>()?;
        let element = Element::from_atomic_number(record[0])
            .ok_or_else(|| self.invalid(format!("unknown atomic number {}", record[0])))?;
        let isotope = NonZeroU16::new(u16::from_le_bytes([record[1], record[2]]));
        let charge = i8::from_le_bytes([record[3]]);
        let radical = match (record[4], record[5]) {
            (0, 0) => None,
            (electrons, spin) => Some(
                AtomRadical::new(electrons, (spin != 0).then_some(spin))
                    .ok_or_else(|| self.invalid("invalid radical"))?,
            ),
        };
        let hydrogens = match (record[6], record[7]) {
            (0, 0) => ImplicitHydrogens::Inferred,
            (1, count) => ImplicitHydrogens::Fixed(count),
            _ => return Err(self.invalid("invalid hydrogen declaration")),
        };
        let atom_map = NonZeroU32::new(u32::from_le_bytes([
            record[8], record[9], record[10], record[11],
        ]));
        Ok(Atom::new(element)
            .with_isotope(isotope)
            .with_formal_charge(charge)
            .with_radical(radical)
            .with_hydrogens(hydrogens)
            .with_atom_map(atom_map))
    }

    fn carrier(&mut self) -> Result<StereoCarrier> {
        match self.u8()? {
            0 => Ok(StereoCarrier::Atom(AtomId::new(self.u32()?))),
            1 => Ok(StereoCarrier::ImplicitHydrogen),
            2 => Ok(StereoCarrier::ImplicitLonePair),
            other => Err(self.invalid(format!("unknown stereo carrier {other}"))),
        }
    }

    fn orientation<T>(&mut self, first: T, second: T) -> Result<Option<T>> {
        match self.u8()? {
            0 => Ok(None),
            1 => Ok(Some(first)),
            2 => Ok(Some(second)),
            other => Err(self.invalid(format!("unknown stereo orientation {other}"))),
        }
    }

    fn stereo_element(&mut self) -> Result<StereoElementKind> {
        match self.u8()? {
            0 => {
                let center = AtomId::new(self.u32()?);
                let carriers = [
                    self.carrier()?,
                    self.carrier()?,
                    self.carrier()?,
                    self.carrier()?,
                ];
                let orientation = self.orientation(
                    TetrahedralOrientation::Clockwise,
                    TetrahedralOrientation::CounterClockwise,
                )?;
                Ok(StereoElementKind::Tetrahedral(TetrahedralStereo {
                    center,
                    carriers,
                    orientation,
                }))
            }
            1 => {
                let bond = BondId::new(self.u32()?);
                let left = AtomId::new(self.u32()?);
                let right = AtomId::new(self.u32()?);
                let left_carrier = self.carrier()?;
                let right_carrier = self.carrier()?;
                let orientation = self.orientation(
                    DoubleBondOrientation::Together,
                    DoubleBondOrientation::Opposite,
                )?;
                Ok(StereoElementKind::DoubleBond(DoubleBondStereo {
                    bond,
                    left,
                    right,
                    left_carrier,
                    right_carrier,
                    orientation,
                }))
            }
            2 => {
                let axis = BondId::new(self.u32()?);
                let carriers = [self.carrier()?, self.carrier()?];
                let orientation = self.orientation(
                    AxisOrientation::Clockwise,
                    AxisOrientation::CounterClockwise,
                )?;
                Ok(StereoElementKind::Axis(AxisStereo {
                    axis,
                    carriers,
                    orientation,
                }))
            }
            other => Err(self.invalid(format!("unknown stereo element kind {other}"))),
        }
    }

    fn atom_ids(&mut self, limit: usize) -> Result<Vec<AtomId>> {
        Ok(self
            .ids(limit, Resource::Atoms)?
            .into_iter()
            .map(AtomId::new)
            .collect())
    }

    fn bond_ids(&mut self, limit: usize) -> Result<Vec<BondId>> {
        Ok(self
            .ids(limit, Resource::Bonds)?
            .into_iter()
            .map(BondId::new)
            .collect())
    }

    fn perception(
        &mut self,
        atom_count: usize,
        bond_count: usize,
        element_count: usize,
    ) -> Result<Perception> {
        let sections = self.u8()?;
        if sections & !0b11_1111 != 0 {
            return Err(self.error(
                DiagnosticKind::Unsupported,
                format!("unknown perception sections {sections:#04x}"),
            ));
        }
        let built = |decoder: &Self, error| {
            decoder
                .invalid("invalid stored perception")
                .with_source(error)
        };
        let mut builder = Perception::builder();
        if sections & 1 != 0 {
            let model = match self.u8()? {
                0 => ValenceModel::RdkitLike,
                other => return Err(self.invalid(format!("unknown valence model {other}"))),
            };
            let validation = match self.u8()? {
                0 => ValenceValidation::Strict,
                1 => ValenceValidation::Permissive,
                other => return Err(self.invalid(format!("unknown valence validation {other}"))),
            };
            let count = self.count(atom_count, Resource::Atoms)?;
            let mut hydrogens = Vec::with_capacity(count);
            for _ in 0..count {
                let atom = AtomId::new(self.u32()?);
                hydrogens.push((atom, self.u8()?));
            }
            builder = builder
                .with_valence(model, validation, hydrogens)
                .map_err(|error| built(self, error))?;
        }
        if sections & 2 != 0 {
            let atom_flags = self.bits(atom_count)?;
            let bond_flags = self.bits(bond_count)?;
            let basis = if self.tag("ring basis")? {
                let model = match self.u8()? {
                    0 => RingBasisModel::FiguerasSssrLike,
                    other => return Err(self.invalid(format!("unknown ring basis model {other}"))),
                };
                let fallback = self.bool()?;
                let ring_count = self.count(bond_count, Resource::Values)?;
                let mut rings = Vec::with_capacity(ring_count.min(CHUNK_BYTES));
                for _ in 0..ring_count {
                    let atoms = self.atom_ids(atom_count)?;
                    let bonds = self.bond_ids(bond_count)?;
                    rings.push(Ring { atoms, bonds });
                }
                let rings = RingSet::from_rings(rings);
                Some(if fallback {
                    RingBasisState::with_depth_first_fallback(model, rings)
                } else {
                    RingBasisState::new(model, rings)
                })
            } else {
                None
            };
            builder = builder.with_rings(RingMembership::from_flags(atom_flags, bond_flags), basis);
        }
        if sections & 4 != 0 {
            let model = match self.u8()? {
                0 => AromaticityModel::RdkitLike,
                1 => AromaticityModel::Mdl,
                other => return Err(self.invalid(format!("unknown aromaticity model {other}"))),
            };
            let atoms = self.atom_ids(atom_count)?;
            let bonds = self.bond_ids(bond_count)?;
            builder = builder
                .with_aromaticity(model, atoms, bonds)
                .map_err(|error| built(self, error))?;
        }
        if sections & 8 != 0 {
            let model = match self.u8()? {
                0 => ConjugationModel::RdkitLike,
                other => return Err(self.invalid(format!("unknown conjugation model {other}"))),
            };
            let atoms = self.atom_ids(atom_count)?;
            let bonds = self.bond_ids(bond_count)?;
            builder = builder
                .with_conjugation(model, atoms, bonds)
                .map_err(|error| built(self, error))?;
        }
        if sections & 16 != 0 {
            let group_count = self.count(bond_count, Resource::Values)?;
            let mut groups = Vec::with_capacity(group_count.min(CHUNK_BYTES));
            for _ in 0..group_count {
                let atoms = self.atom_ids(atom_count)?;
                let bonds = self.bond_ids(bond_count)?;
                groups.push(ResonanceGroup { atoms, bonds });
            }
            builder = builder
                .with_resonance_groups(groups)
                .map_err(|error| built(self, error))?;
        }
        if sections & 32 != 0 {
            let count = self.count(element_count, Resource::Values)?;
            let mut descriptors = Vec::with_capacity(count);
            for _ in 0..count {
                let element = StereoElementId::new(self.u32()?);
                let code = self.u8()?;
                let descriptor = descriptor(code)
                    .ok_or_else(|| self.invalid(format!("unknown CIP descriptor {code}")))?;
                descriptors.push((element, descriptor));
            }
            builder = builder
                .with_cip_descriptors(descriptors)
                .map_err(|error| built(self, error))?;
        }
        Ok(builder.build())
    }

    // Properties.

    fn unit(&mut self) -> Result<Unit> {
        let mut dimension = Dimension::DIMENSIONLESS;
        for base in BASE_DIMENSIONS {
            dimension = dimension.with_exponent(base, self.i32()?);
        }
        let scale = self.f64()?;
        let symbol = self.opt_str()?;
        unit_with_symbol(dimension, scale, symbol.as_deref())
            .map_err(|error| self.invalid("invalid unit").with_source(error))
    }

    fn key(&mut self, previous: Option<&PropertyKey>) -> Result<PropertyKey> {
        let text = self.str()?;
        let key = PropertyKey::new(text)
            .map_err(|error| self.invalid("invalid property key").with_source(error))?;
        if previous.is_some_and(|previous| previous >= &key) {
            return Err(self.invalid("property keys are not strictly ascending"));
        }
        Ok(key)
    }

    fn owner(&mut self) -> Result<OwnerProperties> {
        let count = self.count(self.options.max_items(), Resource::Values)?;
        let mut owner = OwnerProperties::new();
        let mut previous = None;
        for _ in 0..count {
            let key = self.key(previous.as_ref())?;
            let value = match self.u8()? {
                0 => PropertyValue::Bool(self.bool()?),
                1 => PropertyValue::Int(self.i64()?),
                2 => {
                    let value = self.f64()?;
                    let unit = self.unit()?;
                    PropertyValue::Real { value, unit }
                }
                3 => PropertyValue::String(self.str()?),
                other => return Err(self.invalid(format!("unknown property type {other}"))),
            };
            owner
                .insert(key.clone(), value)
                .map_err(|error| self.invalid("invalid property value").with_source(error))?;
            previous = Some(key);
        }
        Ok(owner)
    }

    fn cells<T>(
        &mut self,
        present: &[bool],
        mut read: impl FnMut(&mut Self) -> Result<T>,
    ) -> Result<Vec<Option<T>>> {
        if !present.contains(&true) {
            return Err(self.invalid("a stored column has no values"));
        }
        let mut values = Vec::with_capacity(present.len());
        for &present in present {
            values.push(if present { Some(read(self)?) } else { None });
        }
        Ok(values)
    }

    fn table(&mut self, rows: usize) -> Result<RawPropertyTable> {
        let stored_rows = self.u32()? as usize;
        if stored_rows != rows {
            return Err(self.error(
                DiagnosticKind::Inconsistent,
                format!("a property table has {stored_rows} rows instead of {rows}"),
            ));
        }
        let count = self.count(self.options.max_items(), Resource::Fields)?;
        let mut table = RawPropertyTable::new(rows);
        let mut previous = None;
        for _ in 0..count {
            let key = self.key(previous.as_ref())?;
            let column = match self.u8()? {
                0 => {
                    let present = self.bits(rows)?;
                    PropertyColumn::Bool(self.cells(&present, Self::bool)?)
                }
                1 => {
                    let present = self.bits(rows)?;
                    PropertyColumn::Int(self.cells(&present, Self::i64)?)
                }
                2 => {
                    let unit = self.unit()?;
                    let present = self.bits(rows)?;
                    PropertyColumn::Real {
                        unit,
                        values: self.cells(&present, Self::f64)?,
                    }
                }
                3 => {
                    let present = self.bits(rows)?;
                    PropertyColumn::String(self.cells(&present, Self::str)?)
                }
                other => return Err(self.invalid(format!("unknown column type {other}"))),
            };
            table
                .insert(key.clone(), column)
                .map_err(|error| self.invalid("invalid property column").with_source(error))?;
            previous = Some(key);
        }
        Ok(table)
    }

    // Topology.

    fn topology(&mut self) -> Result<Topology> {
        let max_atoms = self.options.max_atoms();
        let max_bonds = self.options.max_bonds();
        let definition_count = self.count(max_atoms, Resource::Components)?;
        let mut definitions = Vec::with_capacity(definition_count.min(CHUNK_BYTES));
        for _ in 0..definition_count {
            let molecule = self.molecule()?;
            let code = self.u8()?;
            let class = molecule_class(code)
                .ok_or_else(|| self.invalid(format!("unknown molecule class {code}")))?;
            let explicit = self.bool()?;
            definitions.push(RawDefinition {
                molecule,
                class,
                explicit,
            });
        }
        let instance_count = self.count(max_atoms, Resource::Components)?;
        let raw_instances = self.u32_array(instance_count)?;
        let mut instances = Vec::with_capacity(instance_count);
        let (mut atom_total, mut bond_total) = (0usize, 0usize);
        for raw in raw_instances {
            let definition = definitions.get(raw as usize).ok_or_else(|| {
                self.error(
                    DiagnosticKind::InvalidReference,
                    format!("an instance names missing definition {raw}"),
                )
            })?;
            atom_total = atom_total.saturating_add(definition.molecule.atom_count());
            bond_total = bond_total.saturating_add(definition.molecule.bond_count());
            if atom_total > max_atoms {
                return Err(self.limit(Resource::Atoms, max_atoms, atom_total));
            }
            if bond_total > max_bonds {
                return Err(self.limit(Resource::Bonds, max_bonds, bond_total));
            }
            instances.push(MoleculeDefinitionId::new(raw));
        }
        let atom_order = if self.tag("atom order")? {
            let raw = self.u32_array(atom_total.saturating_mul(2))?;
            Some(
                raw.as_chunks::<2>()
                    .0
                    .iter()
                    .map(|pair| {
                        InstanceAtomId::new(MoleculeInstanceId::new(pair[0]), AtomId::new(pair[1]))
                    })
                    .collect(),
            )
        } else {
            None
        };
        // Chains and residues may be empty, so the atom count does not bound
        // them; sites and residue site lists are bounded by the atoms.
        let max_items = self.options.max_items();
        let chain_count = self.count(max_items, Resource::Residues)?;
        let labels = self.strings(chain_count, true)?;
        let authors = self.strings(chain_count, false)?;
        let lists = self.lists(chain_count, max_items, Resource::Residues, ResidueId::new)?;
        let chains = labels
            .into_iter()
            .zip(authors)
            .zip(lists)
            .map(|((label_id, author_id), residues)| RawChain {
                label_id: label_id.unwrap_or_default(),
                author_id,
                residues,
            })
            .collect::<Vec<_>>();
        let residue_count = self.count(max_items, Resource::Residues)?;
        let residues = self.residues(residue_count, atom_total)?;
        let site_count = self.count(atom_total, Resource::Atoms)?;
        let atoms = self.u32_array(site_count.saturating_mul(2))?;
        let type_symbols = self.strings(site_count, false)?;
        let label_asyms = self.strings(site_count, false)?;
        let auth_asyms = self.strings(site_count, false)?;
        let label_atoms = self.strings(site_count, false)?;
        let auth_atoms = self.strings(site_count, false)?;
        let sites = atoms
            .as_chunks::<2>()
            .0
            .iter()
            .zip(type_symbols)
            .zip(label_asyms)
            .zip(auth_asyms)
            .zip(label_atoms)
            .zip(auth_atoms)
            .map(
                |(((((atom, type_symbol), label_asym), auth_asym), label_atom), auth_atom)| {
                    RawAtomSite {
                        atom: InstanceAtomId::new(
                            MoleculeInstanceId::new(atom[0]),
                            AtomId::new(atom[1]),
                        ),
                        metadata: AtomSiteMetadata::default()
                            .with_type_symbol(type_symbol)
                            .with_label_asym_id(label_asym)
                            .with_auth_asym_id(auth_asym)
                            .with_label_atom_id(label_atom)
                            .with_auth_atom_id(auth_atom),
                    }
                },
            )
            .collect::<Vec<_>>();
        let mut properties = TopologyProperties::default();
        *properties.owner_mut() = self.owner()?;
        *properties.molecule_instances_mut() = self.table(instance_count)?;
        *properties.atoms_mut() = self.table(atom_total)?;
        *properties.bonds_mut() = self.table(bond_total)?;
        *properties.chains_mut() = self.table(chain_count)?;
        *properties.residues_mut() = self.table(residue_count)?;
        *properties.atom_sites_mut() = self.table(site_count)?;
        let hierarchy = Hierarchy::from_raw_parts(chains, residues, sites)
            .map_err(|error| self.hierarchy_error(error))?;
        Topology::from_parts(TopologyParts {
            definitions,
            instances,
            atom_order,
            hierarchy,
            properties,
        })
        .map_err(|error| {
            self.rejected(
                DiagnosticKind::Inconsistent,
                "stored topology does not publish",
            )
            .with_source(error)
        })
    }

    fn residues(&mut self, count: usize, site_limit: usize) -> Result<Vec<RawResidue>> {
        let names = self.strings(count, true)?;
        let label_comp_ids = self.strings(count, false)?;
        let author_comp_ids = self.strings(count, false)?;
        let has_label_seq = self.bits(count)?;
        let mut label_seq_ids = Vec::with_capacity(count);
        for present in has_label_seq {
            label_seq_ids.push(if present { Some(self.i32()?) } else { None });
        }
        let author_seq_ids = self.strings(count, false)?;
        let insertion_codes = self.strings(count, false)?;
        let mut codes = Vec::with_capacity(count.min(CHUNK_BYTES));
        self.chunks(count, |bytes| {
            codes.extend_from_slice(bytes);
            Ok(())
        })?;
        let mut classes = Vec::with_capacity(count);
        for code in codes {
            classes.push(
                residue_class(code)
                    .ok_or_else(|| self.invalid(format!("unknown residue class {code}")))?,
            );
        }
        let explicit = self.bits(count)?;
        let sites = self.lists(count, site_limit, Resource::Atoms, AtomSiteId::new)?;
        let mut residues = Vec::with_capacity(count);
        let mut columns = names
            .into_iter()
            .zip(label_comp_ids)
            .zip(author_comp_ids)
            .zip(author_seq_ids)
            .zip(insertion_codes);
        for (index, sites) in sites.into_iter().enumerate() {
            let Some(((((name, label_comp_id), author_comp_id), author_seq_id), insertion_code)) =
                columns.next()
            else {
                break;
            };
            residues.push(RawResidue {
                name: name.unwrap_or_default(),
                label_comp_id,
                author_comp_id,
                label_seq_id: label_seq_ids[index],
                author_seq_id,
                insertion_code,
                class: classes[index],
                explicit: explicit[index],
                atom_sites: sites,
            });
        }
        Ok(residues)
    }

    fn hierarchy_error(&self, error: RawHierarchyError) -> FormatError {
        let kind = match error {
            RawHierarchyError::IdentifierCapacityExceeded(_) => DiagnosticKind::Limit,
            RawHierarchyError::InvalidResidue(_) | RawHierarchyError::InvalidAtomSite(_) => {
                DiagnosticKind::InvalidReference
            }
            _ => DiagnosticKind::Inconsistent,
        };
        self.rejected(kind, format!("stored hierarchy is inconsistent: {error}"))
    }

    // Structures.

    fn points(&mut self, count: usize) -> Result<Vec<Point3>> {
        let mut points = Vec::with_capacity(count);
        self.chunks(count.saturating_mul(24), |bytes| {
            points.extend(bytes.as_chunks::<24>().0.iter().map(|raw| {
                let value = |start: usize| {
                    let mut bytes = [0; 8];
                    bytes.copy_from_slice(&raw[start..start + 8]);
                    f64::from_bits(u64::from_le_bytes(bytes))
                };
                Point3::new(value(0), value(8), value(16))
            }));
            Ok(())
        })?;
        Ok(points)
    }

    fn vectors(&mut self, count: usize) -> Result<Vec<Vector3>> {
        Ok(self
            .points(count)?
            .into_iter()
            .map(|point| Vector3::new(point.x, point.y, point.z))
            .collect())
    }

    fn positions(&mut self, count: usize) -> Result<Positions> {
        let points = self.points(count)?;
        Positions::from_vec(Quantity::new(points, CANONICAL_LENGTH_UNIT))
            .map_err(|error| self.invalid("invalid positions").with_source(error))
    }

    fn cell(&mut self) -> Result<Option<PeriodicCell>> {
        if !self.tag("cell")? {
            return Ok(None);
        }
        let vectors = self.vectors(3)?;
        let axes = self.u8()?;
        if axes & !0b111 != 0 {
            return Err(self.invalid(format!("invalid periodic axes {axes:#04x}")));
        }
        let vectors = [vectors[0], vectors[1], vectors[2]];
        PeriodicCell::new(
            Quantity::new(vectors, CANONICAL_LENGTH_UNIT),
            [axes & 1 != 0, axes & 2 != 0, axes & 4 != 0],
        )
        .map(Some)
        .map_err(|error| self.invalid("invalid periodic cell").with_source(error))
    }

    fn scalars(&mut self, count: usize) -> Result<Option<Vec<Option<f64>>>> {
        if !self.tag("site values")? {
            return Ok(None);
        }
        let present = self.bits(count)?;
        self.cells(&present, Self::f64).map(Some)
    }

    fn structure_properties(&mut self, atoms: usize, bonds: usize) -> Result<StructureProperties> {
        let mut properties = StructureProperties::default();
        *properties.owner_mut() = self.owner()?;
        *properties.atoms_raw_mut() = self.table(atoms)?;
        *properties.bonds_raw_mut() = self.table(bonds)?;
        Ok(properties)
    }

    fn structure_error(
        &self,
        error: impl std::error::Error + Send + Sync + 'static,
    ) -> FormatError {
        self.rejected(DiagnosticKind::InvalidValue, "invalid stored structure")
            .with_source(error)
    }

    fn model(&mut self) -> Result<Model> {
        let topology = Arc::new(self.topology()?);
        let (atoms, bonds) = (topology.atom_count(), topology.bond_count());
        let positions = self.positions(atoms)?;
        let cell = self.cell()?;
        let occupancies = self.scalars(atoms)?;
        let b_factors = self.scalars(atoms)?;
        let properties = self.structure_properties(atoms, bonds)?;
        let mut model =
            Model::new(topology, positions).map_err(|error| self.structure_error(error))?;
        model.set_cell(cell);
        model
            .set_occupancies(occupancies)
            .map_err(|error| self.structure_error(error))?;
        model
            .set_b_factors(b_factors.map(|values| Quantity::new(values, SQUARE_NANOMETER)))
            .map_err(|error| self.structure_error(error))?;
        model
            .set_properties(properties)
            .map_err(|error| self.structure_error(error))?;
        Ok(model)
    }

    fn ensemble(&mut self) -> Result<Ensemble> {
        let topology = Arc::new(self.topology()?);
        let (atoms, bonds) = (topology.atom_count(), topology.bond_count());
        let owner = self.owner()?;
        let temperature = self.opt_f64()?;
        let count = self.count(self.options.max_items(), Resource::Records)?;
        let mut samples = Vec::with_capacity(count.min(CHUNK_BYTES));
        let mut properties = Vec::with_capacity(count.min(CHUNK_BYTES));
        let mut assignments = Vec::with_capacity(count.min(CHUNK_BYTES));
        for _ in 0..count {
            let positions = self.positions(atoms)?;
            let cell = self.cell()?;
            let occupancies = self.scalars(atoms)?;
            let b_factors = self.scalars(atoms)?;
            properties.push(self.structure_properties(atoms, bonds)?);
            let weight = self.f64()?;
            assignments.push(self.u32()? as usize);
            let mut sample = EnsembleSample::new(positions, weight)
                .map_err(|error| self.structure_error(error))?;
            sample.set_cell(cell);
            sample
                .set_occupancies(occupancies)
                .map_err(|error| self.structure_error(error))?;
            sample
                .set_b_factors(b_factors.map(|values| Quantity::new(values, SQUARE_NANOMETER)))
                .map_err(|error| self.structure_error(error))?;
            samples.push(sample);
        }
        let state_count = self.count(self.options.max_items(), Resource::Records)?;
        let mut representatives = Vec::with_capacity(state_count.min(CHUNK_BYTES));
        let mut states = Vec::with_capacity(state_count.min(CHUNK_BYTES));
        for _ in 0..state_count {
            representatives.push(self.u32()? as usize);
            let energy = self.opt_f64()?;
            let mut state = EnsembleState::new();
            state
                .set_energy(energy.map(|energy| Quantity::new(energy, CANONICAL_ENERGY_UNIT)))
                .map_err(|error| self.structure_error(error))?;
            *state.properties_mut() = self.owner()?;
            states.push(state);
        }
        let partition = EnsemblePartition::new(assignments, representatives)
            .map_err(|error| self.structure_error(error))?;
        let mut ensemble = Ensemble::from_samples(topology, samples, partition)
            .map_err(|error| self.structure_error(error))?;
        *ensemble.properties_mut() = owner;
        ensemble
            .set_temperature(
                temperature
                    .map(|temperature| Quantity::new(temperature, CANONICAL_TEMPERATURE_UNIT)),
            )
            .map_err(|error| self.structure_error(error))?;
        for (index, state) in states.into_iter().enumerate() {
            *ensemble
                .state_mut(index)
                .expect("the partition numbers every decoded state") = state;
        }
        for (index, properties) in properties.into_iter().enumerate() {
            ensemble
                .get_mut(index)
                .expect("every decoded sample was added")
                .set_properties(properties)
                .map_err(|error| self.structure_error(error))?;
        }
        Ok(ensemble)
    }

    fn trajectory(&mut self) -> Result<Trajectory> {
        let topology = Arc::new(self.topology()?);
        let (atoms, bonds) = (topology.atom_count(), topology.bond_count());
        let owner = self.owner()?;
        let count = self.count(self.options.max_items(), Resource::Frames)?;
        let mut frames = Vec::with_capacity(count.min(CHUNK_BYTES));
        let mut properties = Vec::with_capacity(count.min(CHUNK_BYTES));
        for _ in 0..count {
            let positions = self.positions(atoms)?;
            let cell = self.cell()?;
            properties.push(self.structure_properties(atoms, bonds)?);
            let velocities = if self.tag("velocities")? {
                let values = self.vectors(atoms)?;
                Some(
                    Velocities::from_vec(Quantity::new(values, CANONICAL_VELOCITY_UNIT))
                        .map_err(|error| self.structure_error(error))?,
                )
            } else {
                None
            };
            let forces = if self.tag("forces")? {
                let values = self.vectors(atoms)?;
                Some(
                    Forces::from_vec(Quantity::new(values, CANONICAL_FORCE_UNIT))
                        .map_err(|error| self.structure_error(error))?,
                )
            } else {
                None
            };
            let time = self.opt_f64()?;
            let step = if self.tag("step")? {
                Some(self.u64()?)
            } else {
                None
            };
            let mut frame = TrajectoryFrame::new(positions);
            frame.set_cell(cell);
            frame
                .set_velocities(velocities)
                .map_err(|error| self.structure_error(error))?;
            frame
                .set_forces(forces)
                .map_err(|error| self.structure_error(error))?;
            frame
                .set_time(time.map(|time| Quantity::new(time, CANONICAL_TIME_UNIT)))
                .map_err(|error| self.structure_error(error))?;
            frame.set_step(step);
            frames.push(frame);
        }
        let mut trajectory = Trajectory::from_items(topology, frames)
            .map_err(|error| self.structure_error(error))?;
        *trajectory.properties_mut() = owner;
        for (index, properties) in properties.into_iter().enumerate() {
            trajectory
                .get_mut(index)
                .expect("every decoded frame was added")
                .set_properties(properties)
                .map_err(|error| self.structure_error(error))?;
        }
        Ok(trajectory)
    }

    fn series(&mut self) -> Result<Series> {
        let topology = Arc::new(self.topology()?);
        let (atoms, bonds) = (topology.atom_count(), topology.bond_count());
        let owner = self.owner()?;
        let axis = match self.u8()? {
            0 => SeriesAxis::Step,
            1 => SeriesAxis::Coordinate(self.unit()?),
            other => return Err(self.invalid(format!("unknown series axis {other}"))),
        };
        let count = self.count(self.options.max_items(), Resource::Records)?;
        let mut items = Vec::with_capacity(count.min(CHUNK_BYTES));
        let mut properties = Vec::with_capacity(count.min(CHUNK_BYTES));
        for _ in 0..count {
            let positions = self.positions(atoms)?;
            let cell = self.cell()?;
            properties.push(self.structure_properties(atoms, bonds)?);
            let item = match axis {
                SeriesAxis::Coordinate(unit) => {
                    SeriesItem::at_coordinate(positions, Quantity::new(self.f64()?, unit))
                        .map_err(|error| self.structure_error(error))?
                }
                SeriesAxis::Step => SeriesItem::at_step(positions, self.u64()?),
            };
            items.push(match cell {
                Some(cell) => item.with_cell(cell),
                None => item,
            });
        }
        let mut series = Series::from_items(topology, axis, items)
            .map_err(|error| self.structure_error(error))?;
        *series.properties_mut() = owner;
        for (index, properties) in properties.into_iter().enumerate() {
            series
                .get_mut(index)
                .expect("every decoded item was added")
                .set_properties(properties)
                .map_err(|error| self.structure_error(error))?;
        }
        Ok(series)
    }
}

fn bond_order(code: u8) -> Option<BondOrder> {
    Some(match code {
        0 => BondOrder::Zero,
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Quadruple,
        5 => BondOrder::Dative,
        _ => return None,
    })
}

fn group_kind(code: u8) -> Option<StereoGroupKind> {
    Some(match code {
        0 => StereoGroupKind::Absolute,
        1 => StereoGroupKind::Relative,
        2 => StereoGroupKind::Racemic,
        3 => StereoGroupKind::And,
        4 => StereoGroupKind::Or,
        _ => return None,
    })
}

fn descriptor(code: u8) -> Option<StereoDescriptor> {
    Some(match code {
        0 => StereoDescriptor::R,
        1 => StereoDescriptor::S,
        2 => StereoDescriptor::LowerR,
        3 => StereoDescriptor::LowerS,
        4 => StereoDescriptor::SeqTrans,
        5 => StereoDescriptor::SeqCis,
        6 => StereoDescriptor::E,
        7 => StereoDescriptor::Z,
        8 => StereoDescriptor::M,
        9 => StereoDescriptor::P,
        10 => StereoDescriptor::LowerM,
        11 => StereoDescriptor::LowerP,
        _ => return None,
    })
}

fn molecule_class(code: u8) -> Option<MoleculeClass> {
    Some(match code {
        0 => MoleculeClass::Protein,
        1 => MoleculeClass::Dna,
        2 => MoleculeClass::Rna,
        3 => MoleculeClass::Carbohydrate,
        4 => MoleculeClass::Water,
        5 => MoleculeClass::Ion,
        6 => MoleculeClass::SmallMolecule,
        7 => MoleculeClass::Other,
        _ => return None,
    })
}

fn residue_class(code: u8) -> Option<ResidueClass> {
    Some(match code {
        0 => ResidueClass::AminoAcid,
        1 => ResidueClass::DnaNucleotide,
        2 => ResidueClass::RnaNucleotide,
        3 => ResidueClass::Carbohydrate,
        4 => ResidueClass::Water,
        5 => ResidueClass::Ion,
        6 => ResidueClass::Other,
        _ => return None,
    })
}
