use super::neighbors::Neighbors;
use super::*;
use crate::{
    core::{
        Atom, AtomId, BondOrder, Element, HydrogenDeclaration, Molecule, MoleculeEditor,
        VanDerWaalsRadiusSource,
    },
    geometry::{PeriodicCell, PeriodicGeometry, Vector3},
    structure::{ModelBuilder, Positions},
    topology::{AtomSiteMetadata, ChainId, MoleculeClass, MoleculeDefinitionId, ResidueClass},
};

// TIP3P oxygen LJ sigma * 2^(1/6)/2, matching OpenMM's exclusion convention.
const WATER_RADIUS: f64 = 0.315_075_240_657_512_4 * 0.5612310241546865;
const SEED_WIDTH: f64 = 3.0;
const DATA: &[u8; 895 * 36] = include_bytes!("data/tip3p.bin");

fn template(index: usize) -> [Point3; 3] {
    std::array::from_fn(|atom| {
        let values: [f64; 3] = std::array::from_fn(|axis| {
            let start = index * 36 + atom * 12 + axis * 4;
            f64::from(i32::from_le_bytes(
                DATA[start..start + 4]
                    .try_into()
                    .expect("fixed asset layout"),
            )) / 10_000.0
        });
        Point3::new(values[0], values[1], values[2])
    })
}
#[derive(Clone, Copy)]
struct Water {
    oxygen: Point3,
    template: usize,
}

pub(super) fn solvate(
    model: &mut Model,
    options: &SolventOptions,
) -> Result<SolvationReport, SolvationError> {
    let cell = *model.cell().ok_or(SolvationError::MissingCell)?;
    if cell.periodic_axes() != [true; 3] {
        return Err(SolvationError::NotFullyPeriodic);
    }
    let geometry = PeriodicGeometry::new(cell)?;
    let shortest = geometry.shortest_translation()?.norm();
    if shortest < WATER_RADIUS {
        return Err(SolvationError::NoSolventSpace);
    }
    let ionic_strength = options.ionic_strength.into_unit(MOLAR)?.into_value();
    if !ionic_strength.is_finite() || ionic_strength < 0.0 {
        return Err(SolvationError::InvalidIonicStrength);
    }
    let ion_separation = options
        .ion_separation
        .into_unit(CANONICAL_LENGTH_UNIT)?
        .into_value();
    if !ion_separation.is_finite() || ion_separation <= 0.0 {
        return Err(SolvationError::InvalidIonSeparation);
    }
    let radii = match &options.solute_radii {
        Some(radii) => {
            if !Arc::ptr_eq(&radii.topology, &model.shared_topology()) {
                return Err(SolvationError::RadiusTopologyMismatch);
            }
            radii.radii.clone()
        }
        None => model
            .atoms()
            .map(|(id, atom)| {
                atom.element
                    .van_der_waals_radius_angstrom(VanDerWaalsRadiusSource::Reference)
                    .map(|r| r * 0.1)
                    .ok_or(SolvationError::MissingElementRadius(id))
            })
            .collect::<Result<Vec<_>, _>>()?,
    };
    let charge = options.net_charge.unwrap_or_else(|| {
        model
            .atoms()
            .map(|(_, atom)| i64::from(atom.formal_charge))
            .sum()
    });
    let origin = center(model.positions().values().value())?;
    let solute: Vec<_> = model
        .positions()
        .values()
        .value()
        .iter()
        .map(|p| Point3::origin() + (*p - origin))
        .collect();
    let cutoffs: Vec<_> = radii.iter().map(|r| r + WATER_RADIUS).collect();
    let max_cutoff = cutoffs.iter().copied().fold(WATER_RADIUS, f64::max);
    let mut solute_search = Neighbors::new(&geometry, max_cutoff)?;
    for &point in &solute {
        solute_search.insert(point)?;
    }
    let waters = pack(
        cell,
        &geometry,
        &solute_search,
        &cutoffs,
        options.max_candidates,
    )?;
    if waters.is_empty() {
        return Err(SolvationError::NoSolventSpace);
    }
    let counterions = if options.neutralize {
        usize::try_from(charge.unsigned_abs()).map_err(|_| SolvationError::CapacityOverflow)?
    } else {
        0
    };
    if counterions > waters.len() {
        return Err(SolvationError::InsufficientIonSites {
            requested: counterions,
            placed: waters.len(),
        });
    }
    let pairs = ((waters.len() - counterions) as f64 * ionic_strength / 55.4 + 0.5).floor();
    if !pairs.is_finite() || pairs >= usize::MAX as f64 {
        return Err(SolvationError::CapacityOverflow);
    }
    let pairs = pairs as usize;
    let positive = pairs
        .checked_add(if charge < 0 { counterions } else { 0 })
        .ok_or(SolvationError::CapacityOverflow)?;
    let negative = pairs
        .checked_add(if charge > 0 { counterions } else { 0 })
        .ok_or(SolvationError::CapacityOverflow)?;
    let ion_count = positive
        .checked_add(negative)
        .ok_or(SolvationError::CapacityOverflow)?;
    if ion_count > waters.len() {
        return Err(SolvationError::InsufficientIonSites {
            requested: ion_count,
            placed: waters.len(),
        });
    }
    if ion_count > 0 && shortest <= ion_separation {
        return Err(SolvationError::InsufficientIonSites {
            requested: ion_count,
            placed: 0,
        });
    }
    let mut ion_sites = Neighbors::new(&geometry, ion_separation)?;
    // Existing monatomic ions participate in separation checks, but are never replaced.
    for (i, (_, atom)) in model.atoms().enumerate() {
        let id = model.topology().atom_ids()[i];
        if model.topology().molecule(id.molecule()).is_ok_and(|m| {
            m.molecule().atom_count() == 1
                && (atom.formal_charge != 0 || m.class() == MoleculeClass::Ion)
        }) {
            ion_sites.insert(solute[i])?;
        }
    }
    let mut order: Vec<_> = (0..waters.len()).collect();
    shuffle(&mut order, options.seed);
    let mut selected = Vec::new();
    selected
        .try_reserve(ion_count)
        .map_err(|_| SolvationError::CapacityOverflow)?;
    for index in order {
        if selected.len() == ion_count {
            break;
        }
        let point = waters[index].oxygen;
        if !ion_sites.any(point, |_, d| d <= ion_separation)? {
            ion_sites.insert(point)?;
            selected.push(index);
        }
    }
    if selected.len() != ion_count {
        return Err(SolvationError::InsufficientIonSites {
            requested: ion_count,
            placed: selected.len(),
        });
    }
    let mut is_ion = vec![false; waters.len()];
    for &index in &selected {
        is_ion[index] = true;
    }
    let report = SolvationReport {
        waters_added: waters.len() - ion_count,
        positive_ions_added: positive,
        negative_ions_added: negative,
        solute_charge: charge,
        charge_basis: if options.net_charge.is_some() {
            SolvationChargeBasis::Override
        } else {
            SolvationChargeBasis::FormalCharges
        },
        seed: options.seed,
        cleared_topology_properties: model
            .topology()
            .properties()
            .iter()
            .map(|(k, _)| k.clone())
            .collect(),
        cleared_model_properties: model.properties().iter().map(|(k, _)| k.clone()).collect(),
    };
    let mut builder = model.to_builder();
    if report.waters_added > 0 {
        let water = definition("O", 0, true);
        let atom_ids: Vec<_> = water.atom_ids().collect();
        let def = builder.add_molecule_definition_owned(water)?;
        builder.set_molecule_class(def, MoleculeClass::Water)?;
        let chain = new_chain(&mut builder, "SOL")?;
        for (index, water) in waters.iter().enumerate().filter(|(i, _)| !is_ion[*i]) {
            let source = template(water.template);
            let mut points = source.map(|p| water.oxygen + (p - source[0]));
            for point in &mut points {
                *point = translated(origin, *point)?;
            }
            append(
                &mut builder,
                def,
                &atom_ids,
                &points,
                chain,
                "HOH",
                &["O", "H1", "H2"],
                index,
                ResidueClass::Water,
            )?;
        }
    }
    for (symbol, charge, indices) in [
        (options.positive_ion.symbol(), 1, &selected[..positive]),
        (options.negative_ion.symbol(), -1, &selected[positive..]),
    ] {
        if indices.is_empty() {
            continue;
        }
        let ion = definition(symbol, charge, false);
        let atom_ids: Vec<_> = ion.atom_ids().collect();
        let def = builder.add_molecule_definition_owned(ion)?;
        builder.set_molecule_class(def, MoleculeClass::Ion)?;
        let chain = new_chain(&mut builder, "ION")?;
        for (index, &site) in indices.iter().enumerate() {
            let point = translated(origin, waters[site].oxygen)?;
            append(
                &mut builder,
                def,
                &atom_ids,
                &[point],
                chain,
                &symbol.to_ascii_uppercase(),
                &[symbol],
                index,
                ResidueClass::Ion,
            )?;
        }
    }
    *model = builder.build()?;
    Ok(report)
}

fn pack(
    cell: PeriodicCell,
    geometry: &PeriodicGeometry,
    solute: &Neighbors<'_>,
    cutoffs: &[f64],
    limit: usize,
) -> Result<Vec<Water>, SolvationError> {
    // AABB of the centered fundamental parallelepiped. Cropping is fractional,
    // so explicit rotated and skewed cells require no restricted vector convention.
    let vectors = cell.vectors().into_value();
    let half = Vector3::new(
        vectors.iter().map(|v| v.x.abs()).sum::<f64>() * 0.5,
        vectors.iter().map(|v| v.y.abs()).sum::<f64>() * 0.5,
        vectors.iter().map(|v| v.z.abs()).sum::<f64>() * 0.5,
    );
    let extents = [half.x, half.y, half.z];
    let mut ranges = [(0_i64, 0_i64); 3];
    let mut count = 895_usize;
    for i in 0..3 {
        let low = (-extents[i] / SEED_WIDTH).floor();
        let high = (extents[i] / SEED_WIDTH).floor();
        if !low.is_finite() || !high.is_finite() || low.abs().max(high.abs()) >= 2.0_f64.powi(48) {
            return Err(SolvationError::CandidateLimit);
        }
        ranges[i] = (low as i64, high as i64);
        count = count
            .checked_mul(
                usize::try_from(ranges[i].1 - ranges[i].0 + 1)
                    .map_err(|_| SolvationError::CandidateLimit)?,
            )
            .ok_or(SolvationError::CandidateLimit)?;
    }
    if count > limit {
        return Err(SolvationError::CandidateLimit);
    }
    let mut waters = Vec::new();
    let mut water_search = Neighbors::new(geometry, WATER_RADIUS)?;
    for x in ranges[0].0..=ranges[0].1 {
        for y in ranges[1].0..=ranges[1].1 {
            for z in ranges[2].0..=ranges[2].1 {
                let shift = Vector3::new(
                    x as f64 * SEED_WIDTH,
                    y as f64 * SEED_WIDTH,
                    z as f64 * SEED_WIDTH,
                );
                for index in 0..895 {
                    let oxygen = template(index)[0] + shift;
                    let fractional = geometry.fractional(oxygen - Point3::origin())?;
                    if fractional.iter().any(|v| *v < -0.5 || *v >= 0.5) {
                        continue;
                    }
                    if solute.any(oxygen, |i, d| d < cutoffs[i])?
                        || water_search.any(oxygen, |_, d| d < WATER_RADIUS)?
                    {
                        continue;
                    }
                    waters
                        .try_reserve(1)
                        .map_err(|_| SolvationError::CapacityOverflow)?;
                    water_search.insert(oxygen)?;
                    waters.push(Water {
                        oxygen,
                        template: index,
                    });
                }
            }
        }
    }
    Ok(waters)
}

fn definition(symbol: &str, charge: i8, water: bool) -> Molecule {
    let mut editor = MoleculeEditor::new();
    let mut atom = Atom::new(Element::from_symbol(symbol).expect("built-in element"));
    atom.formal_charge = charge;
    atom.hydrogens = HydrogenDeclaration::Fixed(0);
    let center = editor.add_atom(atom).expect("small built-in graph");
    if water {
        for _ in 0..2 {
            let mut hydrogen = Atom::new(Element::from_symbol("H").expect("hydrogen"));
            hydrogen.hydrogens = HydrogenDeclaration::Fixed(0);
            let h = editor.add_atom(hydrogen).expect("small built-in graph");
            editor
                .add_bond(center, h, BondOrder::Single)
                .expect("built-in water bond");
        }
    }
    editor.finish().expect("connected built-in species")
}
fn new_chain(builder: &mut ModelBuilder, prefix: &str) -> Result<ChainId, SolvationError> {
    let mut suffix = 1_u64;
    loop {
        let label = format!("{prefix}{suffix}");
        if !builder
            .hierarchy()
            .chains()
            .any(|(_, c)| c.label_id() == label)
        {
            return Ok(builder.hierarchy_mut().add_chain(label, None)?);
        }
        suffix += 1;
    }
}
#[allow(clippy::too_many_arguments)]
fn append(
    builder: &mut ModelBuilder,
    definition: MoleculeDefinitionId,
    atoms: &[AtomId],
    points: &[Point3],
    chain: ChainId,
    residue_name: &str,
    names: &[&str],
    index: usize,
    class: ResidueClass,
) -> Result<(), SolvationError> {
    let positions = Positions::new(Quantity::new(points, CANONICAL_LENGTH_UNIT))?;
    let instance = builder.add_instance(definition, &positions)?;
    let residue = builder.hierarchy_mut().add_residue(
        chain,
        residue_name,
        None,
        Some((index + 1).to_string()),
        None,
    )?;
    for (&atom, &name) in atoms.iter().zip(names) {
        builder.hierarchy_mut().add_atom_site(
            residue,
            InstanceAtomId::new(instance, atom),
            AtomSiteMetadata {
                label_atom_id: Some(name.into()),
                ..Default::default()
            },
        )?;
    }
    builder.set_residue_class(residue, class)?;
    Ok(())
}

// SplitMix64 and unbiased bounded draws; independent of external RNG versions.
fn shuffle<T>(values: &mut [T], mut state: u64) {
    for i in (1..values.len()).rev() {
        let bound = i as u64 + 1;
        let threshold = bound.wrapping_neg() % bound;
        let draw = loop {
            state = state.wrapping_add(0x9e3779b97f4a7c15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
            z ^= z >> 31;
            if z >= threshold {
                break (z % bound) as usize;
            }
        };
        values.swap(i, draw);
    }
}

// Reject offsets that cannot retain the local placement to 1e-8 nm. Checking
// every atom also protects internal water geometry from cancellation on translation.
fn translated(origin: Point3, local: Point3) -> Result<Point3, SolvationError> {
    let delta = local - Point3::origin();
    let placed = origin + delta;
    if !placed.is_finite() || ((placed - origin) - delta).norm() > 1e-8 {
        return Err(PeriodicGeometryError::NumericalFailure.into());
    }
    Ok(placed)
}
