//! Native encoder: one object to a `Write` sink through a bounded staging
//! buffer.

use std::io::Write;

use super::crc::Crc32;
use super::{
    NativeContent, NativeKind, NativeSource, NativeWriteOptions, NATIVE_FORMAT_VERSION,
    NATIVE_MAGIC,
};
use crate::diagnostics::{DiagnosticKind, Format, FormatError};
use crate::geometry::{PeriodicCell, Point3, Vector3};
use crate::molecule::{
    AtomRadical, AxisOrientation, BondOrder, DoubleBondOrientation, ImplicitHydrogens, Molecule,
    StereoCarrier, StereoDescriptor, StereoElementKind, StereoGroupKind, TetrahedralOrientation,
};
use crate::perception::{
    aromaticity::AromaticityModel, conjugation::ConjugationModel, rings::RingBasisModel,
    valence::ValenceModel, valence::ValenceValidation, Perception,
};
use crate::properties::{
    OwnerProperties, PropertyColumn, PropertyRow, PropertyTable, PropertyValue, StructureProperties,
};
use crate::structure::{Ensemble, Model, Trajectory};
use crate::topology::{MoleculeClass, ResidueClass, Topology};
use crate::units::{BaseDimension, Unit};

/// The staging buffer size: output reaches the sink in chunks of this size.
const STAGE_BYTES: usize = 64 * 1024;

pub(super) const END_MARKER: u32 = 0x444E_454B;
pub(super) const ABSENT_STRING: u32 = u32::MAX;

pub(crate) fn write_native_to<'a>(
    writer: &mut impl Write,
    source: impl Into<NativeSource<'a>>,
    options: NativeWriteOptions,
) -> Result {
    let source = source.into();
    let mut encoder = Encoder {
        sink: writer,
        buffer: Vec::with_capacity(STAGE_BYTES),
        crc: Crc32::new(),
        perception: options.perception(),
    };
    encoder.header(source.kind())?;
    match source.content {
        NativeContent::Molecule(molecule) => encoder.molecule(molecule)?,
        NativeContent::Topology(topology) => encoder.topology(topology)?,
        NativeContent::Model(model) => encoder.model(model)?,
        NativeContent::Ensemble(ensemble) => encoder.ensemble(ensemble)?,
        NativeContent::Trajectory(trajectory) => encoder.trajectory(trajectory)?,
    }
    encoder.finish()
}

pub(crate) fn write_native<'a>(
    source: impl Into<NativeSource<'a>>,
    options: NativeWriteOptions,
) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    write_native_to(&mut bytes, source, options)?;
    Ok(bytes)
}

struct Encoder<'w, W: Write> {
    sink: &'w mut W,
    buffer: Vec<u8>,
    crc: Crc32,
    perception: bool,
}

fn unrepresentable(message: impl Into<String>) -> FormatError {
    FormatError::write(Format::Native, DiagnosticKind::Unsupported, message)
}

type Result<T = ()> = std::result::Result<T, FormatError>;

impl<W: Write> Encoder<'_, W> {
    fn flush_stage(&mut self) -> Result {
        self.crc.update(&self.buffer);
        self.sink
            .write_all(&self.buffer)
            .map_err(|error| FormatError::sink(Format::Native, error))?;
        self.buffer.clear();
        Ok(())
    }

    fn bytes(&mut self, bytes: &[u8]) -> Result {
        if self.buffer.len() + bytes.len() > STAGE_BYTES {
            self.flush_stage()?;
        }
        if bytes.len() >= STAGE_BYTES {
            self.crc.update(bytes);
            return self
                .sink
                .write_all(bytes)
                .map_err(|error| FormatError::sink(Format::Native, error));
        }
        self.buffer.extend_from_slice(bytes);
        Ok(())
    }

    fn u8(&mut self, value: u8) -> Result {
        self.bytes(&[value])
    }

    fn u16(&mut self, value: u16) -> Result {
        self.bytes(&value.to_le_bytes())
    }

    fn u32(&mut self, value: u32) -> Result {
        self.bytes(&value.to_le_bytes())
    }

    fn u64(&mut self, value: u64) -> Result {
        self.bytes(&value.to_le_bytes())
    }

    fn i32(&mut self, value: i32) -> Result {
        self.bytes(&value.to_le_bytes())
    }

    fn i64(&mut self, value: i64) -> Result {
        self.bytes(&value.to_le_bytes())
    }

    fn f64(&mut self, value: f64) -> Result {
        self.bytes(&value.to_bits().to_le_bytes())
    }

    fn bool(&mut self, value: bool) -> Result {
        self.u8(u8::from(value))
    }

    fn count(&mut self, count: usize) -> Result {
        let count = u32::try_from(count)
            .map_err(|_| unrepresentable(format!("{count} entries exceed a native count")))?;
        self.u32(count)
    }

    fn str(&mut self, value: &str) -> Result {
        if u32::try_from(value.len()).map_or(true, |len| len == ABSENT_STRING) {
            return Err(unrepresentable("a string exceeds the native length range"));
        }
        self.count(value.len())?;
        self.bytes(value.as_bytes())
    }

    fn opt_str(&mut self, value: Option<&str>) -> Result {
        match value {
            Some(value) => self.str(value),
            None => self.u32(ABSENT_STRING),
        }
    }

    fn bits(&mut self, flags: impl ExactSizeIterator<Item = bool>) -> Result {
        let mut byte = 0u8;
        let mut filled = 0;
        for flag in flags {
            byte |= u8::from(flag) << filled;
            filled += 1;
            if filled == 8 {
                self.u8(byte)?;
                byte = 0;
                filled = 0;
            }
        }
        if filled > 0 {
            self.u8(byte)?;
        }
        Ok(())
    }

    /// A string column: every length (absent strings as [`ABSENT_STRING`]),
    /// then every present string's bytes.
    fn strings<'s, I: Iterator<Item = Option<&'s str>>>(
        &mut self,
        values: impl Fn() -> I,
    ) -> Result {
        for value in values() {
            match value {
                Some(value) => {
                    if u32::try_from(value.len()).map_or(true, |len| len == ABSENT_STRING) {
                        return Err(unrepresentable("a string exceeds the native length range"));
                    }
                    self.count(value.len())?;
                }
                None => self.u32(ABSENT_STRING)?,
            }
        }
        for value in values().flatten() {
            self.bytes(value.as_bytes())?;
        }
        Ok(())
    }

    /// A column of ID lists: every length, then every list's IDs.
    fn lists<L: ExactSizeIterator<Item = u32>, I: Iterator<Item = L>>(
        &mut self,
        lists: impl Fn() -> I,
    ) -> Result {
        for list in lists() {
            self.count(list.len())?;
        }
        for id in lists().flatten() {
            self.u32(id)?;
        }
        Ok(())
    }

    fn ids(&mut self, ids: impl ExactSizeIterator<Item = u32>) -> Result {
        self.count(ids.len())?;
        for id in ids {
            self.u32(id)?;
        }
        Ok(())
    }

    fn header(&mut self, kind: NativeKind) -> Result {
        self.bytes(&NATIVE_MAGIC)?;
        self.u32(NATIVE_FORMAT_VERSION)?;
        self.u8(kind.code())?;
        self.u8(u8::from(self.perception))?;
        self.u16(0)
    }

    fn finish(mut self) -> Result {
        self.u32(END_MARKER)?;
        self.flush_stage()?;
        let crc = self.crc.finish();
        self.sink
            .write_all(&crc.to_le_bytes())
            .map_err(|error| FormatError::sink(Format::Native, error))
    }

    // Molecules.

    fn molecule(&mut self, molecule: &Molecule) -> Result {
        self.count(molecule.atom_count())?;
        for (_, atom) in molecule.atoms() {
            self.u8(atom.element().atomic_number())?;
            self.u16(atom.isotope().map_or(0, |isotope| isotope.get()))?;
            self.bytes(&atom.formal_charge().to_le_bytes())?;
            let (electrons, spin) = atom.radical().map_or((0, 0), radical_codes);
            self.u8(electrons)?;
            self.u8(spin)?;
            let (kind, count) = match atom.hydrogens() {
                ImplicitHydrogens::Inferred => (0, 0),
                ImplicitHydrogens::Fixed(count) => (1, count),
            };
            self.u8(kind)?;
            self.u8(count)?;
            self.u32(atom.atom_map().map_or(0, |map| map.get()))?;
        }
        self.count(molecule.bond_count())?;
        for (_, bond) in molecule.bonds() {
            self.u32(bond.a().raw())?;
            self.u32(bond.b().raw())?;
            self.u8(bond_order_code(bond.order))?;
        }
        let adjacency = &molecule.graph().adjacency;
        if adjacency
            .iter()
            .all(|bonds| bonds.windows(2).all(|pair| pair[0] < pair[1]))
        {
            self.u8(0)?;
        } else {
            self.u8(1)?;
            for bonds in adjacency {
                for bond in bonds {
                    self.u32(bond.raw())?;
                }
            }
        }
        self.count(molecule.stereo_elements().count())?;
        for (_, element) in molecule.stereo_elements() {
            self.stereo_element(&element.kind)?;
        }
        self.count(molecule.stereo_groups().count())?;
        for (_, group) in molecule.stereo_groups() {
            self.u8(group_kind_code(group.kind))?;
            self.ids(group.members.iter().map(|member| member.raw()))?;
        }
        let properties = molecule.properties();
        self.owner(properties.owner())?;
        self.table(properties.atoms())?;
        self.table(properties.bonds())?;
        if self.perception {
            self.perception(molecule.perception())?;
        }
        Ok(())
    }

    fn carrier(&mut self, carrier: StereoCarrier) -> Result {
        match carrier {
            StereoCarrier::Atom(atom) => {
                self.u8(0)?;
                self.u32(atom.raw())
            }
            StereoCarrier::ImplicitHydrogen => self.u8(1),
            StereoCarrier::ImplicitLonePair => self.u8(2),
        }
    }

    fn stereo_element(&mut self, kind: &StereoElementKind) -> Result {
        match kind {
            StereoElementKind::Tetrahedral(stereo) => {
                self.u8(0)?;
                self.u32(stereo.center.raw())?;
                for carrier in stereo.carriers {
                    self.carrier(carrier)?;
                }
                self.u8(match stereo.orientation {
                    None => 0,
                    Some(TetrahedralOrientation::Clockwise) => 1,
                    Some(TetrahedralOrientation::CounterClockwise) => 2,
                })
            }
            StereoElementKind::DoubleBond(stereo) => {
                self.u8(1)?;
                self.u32(stereo.bond.raw())?;
                self.u32(stereo.left.raw())?;
                self.u32(stereo.right.raw())?;
                self.carrier(stereo.left_carrier)?;
                self.carrier(stereo.right_carrier)?;
                self.u8(match stereo.orientation {
                    None => 0,
                    Some(DoubleBondOrientation::Together) => 1,
                    Some(DoubleBondOrientation::Opposite) => 2,
                })
            }
            StereoElementKind::Axis(stereo) => {
                self.u8(2)?;
                self.u32(stereo.axis.raw())?;
                for carrier in stereo.carriers {
                    self.carrier(carrier)?;
                }
                self.u8(match stereo.orientation {
                    None => 0,
                    Some(AxisOrientation::Clockwise) => 1,
                    Some(AxisOrientation::CounterClockwise) => 2,
                })
            }
        }
    }

    fn perception(&mut self, perception: &Perception) -> Result {
        let sections = [
            perception.valence().is_some(),
            perception.rings().is_some(),
            perception.aromaticity().is_some(),
            perception.conjugation().is_some(),
            perception.resonance().is_some(),
            perception.cip().is_some(),
        ];
        self.u8(sections
            .iter()
            .enumerate()
            .fold(0, |mask, (bit, present)| mask | (u8::from(*present) << bit)))?;
        if let Some(valence) = perception.valence() {
            self.u8(match valence.model() {
                ValenceModel::RdkitLike => 0,
            })?;
            self.u8(match valence.validation() {
                ValenceValidation::Strict => 0,
                ValenceValidation::Permissive => 1,
            })?;
            self.count(valence.inferred_hydrogens().count())?;
            for (atom, count) in valence.inferred_hydrogens() {
                self.u32(atom.raw())?;
                self.u8(count)?;
            }
        }
        if let Some(rings) = perception.rings() {
            let membership = rings.membership();
            self.bits(membership.atom_flags().iter().copied())?;
            self.bits(membership.bond_flags().iter().copied())?;
            match rings.basis() {
                None => self.u8(0)?,
                Some(basis) => {
                    self.u8(1)?;
                    self.u8(match basis.model() {
                        RingBasisModel::FiguerasSssrLike => 0,
                    })?;
                    self.bool(basis.used_depth_first_fallback())?;
                    let rings = basis.ring_set().rings();
                    self.count(rings.len())?;
                    for ring in rings {
                        self.ids(ring.atoms.iter().map(|atom| atom.raw()))?;
                        self.ids(ring.bonds.iter().map(|bond| bond.raw()))?;
                    }
                }
            }
        }
        if let Some(aromaticity) = perception.aromaticity() {
            self.u8(match aromaticity.model() {
                AromaticityModel::RdkitLike => 0,
                AromaticityModel::Mdl => 1,
            })?;
            self.ids(aromaticity.atoms().map(|atom| atom.raw()))?;
            self.ids(aromaticity.bonds().map(|bond| bond.raw()))?;
        }
        if let Some(conjugation) = perception.conjugation() {
            self.u8(match conjugation.model() {
                ConjugationModel::RdkitLike => 0,
            })?;
            self.ids(conjugation.atoms().map(|atom| atom.raw()))?;
            self.ids(conjugation.bonds().map(|bond| bond.raw()))?;
        }
        if let Some(resonance) = perception.resonance() {
            self.count(resonance.groups().len())?;
            for group in resonance.groups() {
                self.ids(group.atoms.iter().map(|atom| atom.raw()))?;
                self.ids(group.bonds.iter().map(|bond| bond.raw()))?;
            }
        }
        if let Some(cip) = perception.cip() {
            self.count(cip.cip_descriptors().count())?;
            for (element, descriptor) in cip.cip_descriptors() {
                self.u32(element.raw())?;
                self.u8(descriptor_code(descriptor))?;
            }
        }
        Ok(())
    }

    // Properties.

    fn unit(&mut self, unit: Unit) -> Result {
        let dimension = unit.dimension();
        for base in BASE_DIMENSIONS {
            self.i32(dimension.exponent(base))?;
        }
        self.f64(unit.scale())?;
        self.opt_str(unit.symbol())
    }

    fn owner(&mut self, owner: &OwnerProperties) -> Result {
        self.count(owner.len())?;
        for (key, value) in owner.iter() {
            self.str(key.as_str())?;
            match value {
                PropertyValue::Bool(value) => {
                    self.u8(0)?;
                    self.bool(*value)?;
                }
                PropertyValue::Int(value) => {
                    self.u8(1)?;
                    self.i64(*value)?;
                }
                PropertyValue::Real { value, unit } => {
                    self.u8(2)?;
                    self.f64(*value)?;
                    self.unit(*unit)?;
                }
                PropertyValue::String(value) => {
                    self.u8(3)?;
                    self.str(value)?;
                }
            }
        }
        Ok(())
    }

    fn table<R: PropertyRow>(&mut self, table: &PropertyTable<R>) -> Result {
        self.count(table.len())?;
        self.count(table.iter().len())?;
        for (key, column) in table.iter() {
            self.str(key.as_str())?;
            match column {
                PropertyColumn::Bool(values) => {
                    self.u8(0)?;
                    self.bits(values.iter().map(Option::is_some))?;
                    for value in values.iter().flatten() {
                        self.bool(*value)?;
                    }
                }
                PropertyColumn::Int(values) => {
                    self.u8(1)?;
                    self.bits(values.iter().map(Option::is_some))?;
                    for value in values.iter().flatten() {
                        self.i64(*value)?;
                    }
                }
                PropertyColumn::Real { unit, values } => {
                    self.u8(2)?;
                    self.unit(*unit)?;
                    self.bits(values.iter().map(Option::is_some))?;
                    for value in values.iter().flatten() {
                        self.f64(*value)?;
                    }
                }
                PropertyColumn::String(values) => {
                    self.u8(3)?;
                    self.bits(values.iter().map(Option::is_some))?;
                    for value in values.iter().flatten() {
                        self.str(value)?;
                    }
                }
            }
        }
        Ok(())
    }

    // Topology.

    fn topology(&mut self, topology: &Topology) -> Result {
        self.count(topology.definition_count())?;
        for definition in topology.definitions() {
            self.molecule(definition.molecule())?;
            self.u8(molecule_class_code(definition.class()))?;
            self.bool(topology.definition_class_is_explicit(definition.id()))?;
        }
        self.count(topology.instance_count())?;
        for instance in topology.instances() {
            self.u32(instance.definition_id().raw())?;
        }
        let instance_order = topology
            .instances()
            .flat_map(|instance| {
                let id = instance.id();
                instance
                    .molecule()
                    .atom_ids()
                    .map(move |atom| crate::topology::InstanceAtomId::new(id, atom))
            })
            .eq(topology.atom_ids().iter().copied());
        if instance_order {
            self.u8(0)?;
        } else {
            self.u8(1)?;
            for atom in topology.atom_ids() {
                self.u32(atom.instance().raw())?;
                self.u32(atom.atom().raw())?;
            }
        }
        let hierarchy = topology.hierarchy();
        let chains = || hierarchy.chains().map(|(_, chain)| chain);
        self.count(hierarchy.chains().len())?;
        self.strings(|| chains().map(|chain| Some(chain.label_id())))?;
        self.strings(|| chains().map(|chain| chain.author_id()))?;
        self.lists(|| chains().map(|chain| chain.residues().iter().map(|id| id.raw())))?;
        let residues = || hierarchy.residues().map(|(_, residue)| residue);
        self.count(hierarchy.residues().len())?;
        self.strings(|| residues().map(|residue| Some(residue.name())))?;
        self.strings(|| residues().map(|residue| residue.label_comp_id()))?;
        self.strings(|| residues().map(|residue| residue.author_comp_id()))?;
        self.bits(residues().map(|residue| residue.label_seq_id().is_some()))?;
        for sequence in residues().filter_map(|residue| residue.label_seq_id()) {
            self.i32(sequence)?;
        }
        self.strings(|| residues().map(|residue| residue.author_seq_id()))?;
        self.strings(|| residues().map(|residue| residue.insertion_code()))?;
        for residue in residues() {
            self.u8(residue_class_code(residue.class()))?;
        }
        self.bits(residues().map(|residue| residue.explicit_class.is_some()))?;
        self.lists(|| residues().map(|residue| residue.atom_sites().iter().map(|id| id.raw())))?;
        let sites = || hierarchy.atom_sites().map(|(_, site)| site);
        self.count(hierarchy.atom_sites().len())?;
        for site in sites() {
            self.u32(site.atom().instance().raw())?;
            self.u32(site.atom().atom().raw())?;
        }
        self.strings(|| sites().map(|site| site.metadata().type_symbol()))?;
        self.strings(|| sites().map(|site| site.metadata().label_asym_id()))?;
        self.strings(|| sites().map(|site| site.metadata().auth_asym_id()))?;
        self.strings(|| sites().map(|site| site.metadata().label_atom_id()))?;
        self.strings(|| sites().map(|site| site.metadata().auth_atom_id()))?;
        let properties = topology.properties();
        self.owner(properties.owner())?;
        self.table(properties.molecule_instances())?;
        self.table(properties.atoms())?;
        self.table(properties.bonds())?;
        self.table(properties.chains())?;
        self.table(properties.residues())?;
        self.table(properties.atom_sites())
    }

    // Structures.

    fn points(&mut self, points: &[Point3]) -> Result {
        for point in points {
            self.f64(point.x)?;
            self.f64(point.y)?;
            self.f64(point.z)?;
        }
        Ok(())
    }

    fn vectors(&mut self, vectors: &[Vector3]) -> Result {
        for vector in vectors {
            self.f64(vector.x)?;
            self.f64(vector.y)?;
            self.f64(vector.z)?;
        }
        Ok(())
    }

    fn cell(&mut self, cell: Option<&PeriodicCell>) -> Result {
        let Some(cell) = cell else {
            return self.u8(0);
        };
        self.u8(1)?;
        self.vectors(cell.vectors().value())?;
        let axes = cell.periodic_axes();
        self.u8(u8::from(axes[0]) | (u8::from(axes[1]) << 1) | (u8::from(axes[2]) << 2))
    }

    fn scalars(&mut self, values: Option<&[Option<f64>]>) -> Result {
        let Some(values) = values else {
            return self.u8(0);
        };
        self.u8(1)?;
        self.bits(values.iter().map(Option::is_some))?;
        for value in values.iter().flatten() {
            self.f64(*value)?;
        }
        Ok(())
    }

    fn structure_properties(&mut self, properties: &StructureProperties) -> Result {
        self.owner(properties.owner())?;
        self.table(properties.atoms())?;
        self.table(properties.bonds())
    }

    fn model(&mut self, model: &Model) -> Result {
        self.topology(model.topology())?;
        self.points(model.positions().values().value())?;
        self.cell(model.cell())?;
        self.scalars(model.occupancies())?;
        self.scalars(model.b_factors().map(|values| *values.value()))?;
        self.structure_properties(model.properties())
    }

    fn ensemble(&mut self, ensemble: &Ensemble) -> Result {
        self.topology(ensemble.topology())?;
        self.owner(ensemble.properties())?;
        self.count(ensemble.len())?;
        for member in ensemble.iter() {
            let member = member.payload();
            self.points(member.positions().values().value())?;
            self.cell(member.cell())?;
            self.scalars(member.occupancies())?;
            self.scalars(member.b_factors().map(|values| *values.value()))?;
            self.structure_properties(member.properties())?;
            self.f64(member.weight())?;
            self.opt_f64(member.energy().map(|energy| *energy.value()))?;
        }
        Ok(())
    }

    fn trajectory(&mut self, trajectory: &Trajectory) -> Result {
        self.topology(trajectory.topology())?;
        self.owner(trajectory.properties())?;
        self.count(trajectory.len())?;
        for frame in trajectory.iter() {
            let frame = frame.payload();
            self.points(frame.positions().values().value())?;
            self.cell(frame.cell())?;
            self.structure_properties(frame.properties())?;
            match frame.velocities() {
                None => self.u8(0)?,
                Some(velocities) => {
                    self.u8(1)?;
                    self.vectors(velocities.values().value())?;
                }
            }
            match frame.forces() {
                None => self.u8(0)?,
                Some(forces) => {
                    self.u8(1)?;
                    self.vectors(forces.values().value())?;
                }
            }
            self.opt_f64(frame.time().map(|time| *time.value()))?;
            match frame.step() {
                None => self.u8(0)?,
                Some(step) => {
                    self.u8(1)?;
                    self.u64(step)?;
                }
            }
        }
        Ok(())
    }

    fn opt_f64(&mut self, value: Option<f64>) -> Result {
        match value {
            None => self.u8(0),
            Some(value) => {
                self.u8(1)?;
                self.f64(value)
            }
        }
    }
}

pub(super) const BASE_DIMENSIONS: [BaseDimension; 7] = [
    BaseDimension::Length,
    BaseDimension::Mass,
    BaseDimension::Time,
    BaseDimension::Temperature,
    BaseDimension::Amount,
    BaseDimension::Charge,
    BaseDimension::Angle,
];

fn radical_codes(radical: AtomRadical) -> (u8, u8) {
    (
        radical.electron_count(),
        radical.spin_multiplicity().unwrap_or(0),
    )
}

fn bond_order_code(order: BondOrder) -> u8 {
    match order {
        BondOrder::Zero => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Dative => 5,
    }
}

fn group_kind_code(kind: StereoGroupKind) -> u8 {
    match kind {
        StereoGroupKind::Absolute => 0,
        StereoGroupKind::Relative => 1,
        StereoGroupKind::Racemic => 2,
        StereoGroupKind::And => 3,
        StereoGroupKind::Or => 4,
    }
}

fn descriptor_code(descriptor: StereoDescriptor) -> u8 {
    match descriptor {
        StereoDescriptor::R => 0,
        StereoDescriptor::S => 1,
        StereoDescriptor::LowerR => 2,
        StereoDescriptor::LowerS => 3,
        StereoDescriptor::SeqTrans => 4,
        StereoDescriptor::SeqCis => 5,
        StereoDescriptor::E => 6,
        StereoDescriptor::Z => 7,
        StereoDescriptor::M => 8,
        StereoDescriptor::P => 9,
        StereoDescriptor::LowerM => 10,
        StereoDescriptor::LowerP => 11,
    }
}

fn molecule_class_code(class: MoleculeClass) -> u8 {
    match class {
        MoleculeClass::Protein => 0,
        MoleculeClass::Dna => 1,
        MoleculeClass::Rna => 2,
        MoleculeClass::Carbohydrate => 3,
        MoleculeClass::Water => 4,
        MoleculeClass::Ion => 5,
        MoleculeClass::SmallMolecule => 6,
        MoleculeClass::Other => 7,
    }
}

fn residue_class_code(class: ResidueClass) -> u8 {
    match class {
        ResidueClass::AminoAcid => 0,
        ResidueClass::DnaNucleotide => 1,
        ResidueClass::RnaNucleotide => 2,
        ResidueClass::Carbohydrate => 3,
        ResidueClass::Water => 4,
        ResidueClass::Ion => 5,
        ResidueClass::Other => 6,
    }
}
