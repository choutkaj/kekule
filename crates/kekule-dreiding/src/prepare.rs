use std::collections::{BTreeMap, BTreeSet, HashMap};

use dreid_forge::{
    forge, AnglePotential, AnglePotentialType, Atom as ForgeAtom, Bond as ForgeBond,
    BondOrder as ForgeBondOrder, BondPotential, BondPotentialType, ChargeMethod, ForgeConfig,
    ForgedSystem, HBondPotential, InversionPotential, QeqConfig, System, TorsionPotential,
    VdwPairPotential, VdwPotentialType,
};
use kekule::core::BondOrder;
use kekule::structure::ModelView;
use kekule::topology::{
    InstanceAtomId, InstanceBondId, MoleculeInstanceId, Topology, TopologyIdentity,
};
use kekule::units::{Quantity, ELEMENTARY_CHARGE};

use crate::DreidingPrepareError;

pub(crate) const KCAL_TO_KJ: f64 = 4.184;
pub(crate) const COULOMB_KJ_ANGSTROM_PER_MOL_E2: f64 = 1_389.354_576_443_82;

/// Scope used for fixed QEq charge preparation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum QeqGrouping {
    WholeTopology,
    MoleculeInstances,
    ConnectedComponents,
}

/// Explicit DREIDING preparation policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DreidingPrepareOptions {
    pub qeq_grouping: QeqGrouping,
}

impl Default for DreidingPrepareOptions {
    fn default() -> Self {
        Self {
            qeq_grouping: QeqGrouping::MoleculeInstances,
        }
    }
}

/// A prepared, fixed-charge DREIDING potential for one exact topology.
///
/// The reference view is used explicitly for the geometry-dependent QEq
/// preparation performed by `dreid-forge`. Evaluation can then consume any
/// nonperiodic model, ensemble member, or trajectory frame with the same exact
/// topology. Periodic cells are unsupported during preparation and evaluation.
#[derive(Debug, Clone)]
pub struct DreidingPotential {
    pub(crate) topology: TopologyIdentity,
    pub(crate) qeq_grouping: QeqGrouping,
    pub(crate) atom_ids: Vec<InstanceAtomId>,
    pub(crate) atom_indexes: BTreeMap<InstanceAtomId, usize>,
    pub(crate) atom_types: Vec<String>,
    pub(crate) partial_charges: Vec<f64>,
    pub(crate) bonds: Vec<BondTerm>,
    pub(crate) angles: Vec<AngleTerm>,
    pub(crate) torsions: Vec<TorsionTerm>,
    pub(crate) inversions: Vec<InversionTerm>,
    pub(crate) nonbonded: Vec<NonbondedTerm>,
    pub(crate) hydrogen_bonds: Vec<HydrogenBondTerm>,
}

impl DreidingPotential {
    /// Prepares standard DREIDING parameters and fixed QEq charges.
    ///
    /// The reference view must be nonperiodic because the present adapter has
    /// no complete minimum-image or unwrapped-coordinate preparation policy.
    pub fn prepare(
        topology: &Topology,
        reference: ModelView<'_>,
        options: DreidingPrepareOptions,
    ) -> Result<Self, DreidingPrepareError> {
        if !topology.same_identity(reference.topology()) {
            return Err(DreidingPrepareError::TopologyIdentityMismatch);
        }
        if reference.cell().is_some() {
            return Err(DreidingPrepareError::UnsupportedPeriodicCell);
        }
        let prepared = PreparedInput::new(topology, reference, options.qeq_grouping)?;
        let total_charge = prepared
            .groups
            .iter()
            .map(|group| group.formal_charge)
            .sum();
        let whole = forge_system(&prepared.whole, total_charge, None)?;
        let whole_types = per_atom_types(&whole)?;

        let mut partial_charges = vec![0.0; topology.atom_count()];
        for group in &prepared.groups {
            let forged = forge_system(&group.system, group.formal_charge, group.instance)?;
            let local_types = per_atom_types(&forged)?;
            for (local, &global) in group.global_atoms.iter().enumerate() {
                let whole_type = &whole_types[global];
                let group_type = &local_types[local];
                if whole_type != group_type {
                    return Err(DreidingPrepareError::AtomTypeMismatch {
                        molecule: group.instance,
                        atom: topology.atom_ids()[global],
                        whole_model: whole_type.clone(),
                        component_model: group_type.clone(),
                    });
                }
                partial_charges[global] = forged.atom_properties[local].charge;
            }
        }
        require_finite_slice("partial charge", &partial_charges)?;

        let adjacency = adjacency(topology);
        let exclusions = nonbonded_exclusions(&adjacency);
        let bonds = prepare_bonds(&whole)?;
        let angles = prepare_angles(&whole)?;
        let torsions = prepare_torsions(&whole)?;
        let inversions = prepare_inversions(&whole)?;
        let nonbonded = prepare_nonbonded(topology, &whole, &partial_charges, &exclusions)?;
        let hydrogen_bonds = prepare_hydrogen_bonds(&whole, &adjacency, &exclusions)?;

        Ok(Self {
            topology: topology.identity(),
            qeq_grouping: options.qeq_grouping,
            atom_ids: topology.atom_ids().to_vec(),
            atom_indexes: topology
                .atom_ids()
                .iter()
                .enumerate()
                .map(|(index, atom)| (*atom, index))
                .collect(),
            atom_types: whole_types,
            partial_charges,
            bonds,
            angles,
            torsions,
            inversions,
            nonbonded,
            hydrogen_bonds,
        })
    }

    pub const fn qeq_grouping(&self) -> QeqGrouping {
        self.qeq_grouping
    }

    /// Returns the DREIDING type assigned to a model atom.
    pub fn atom_type(&self, atom: InstanceAtomId) -> Option<&str> {
        self.atom_types
            .get(*self.atom_indexes.get(&atom)?)
            .map(String::as_str)
    }

    /// Returns the fixed QEq partial charge assigned to a model atom.
    pub fn partial_charge(&self, atom: InstanceAtomId) -> Option<Quantity<f64>> {
        self.partial_charges
            .get(*self.atom_indexes.get(&atom)?)
            .copied()
            .map(|charge| Quantity::new(charge, ELEMENTARY_CHARGE))
    }
}

#[derive(Debug)]
struct PreparedInput {
    whole: System,
    groups: Vec<PreparedChargeGroup>,
}

impl PreparedInput {
    fn new(
        topology: &Topology,
        reference: ModelView<'_>,
        grouping: QeqGrouping,
    ) -> Result<Self, DreidingPrepareError> {
        validate_atoms(topology)?;
        validate_bonds(topology)?;

        let whole_atoms = topology.atom_ids().to_vec();
        let whole = system_from_selection(topology, reference, &whole_atoms)?;
        let selections = match grouping {
            QeqGrouping::WholeTopology => vec![(None, whole_atoms)],
            QeqGrouping::MoleculeInstances => topology
                .instances()
                .map(|(id, _)| {
                    let atoms = topology
                        .graph_for_instance(id)
                        .expect("validated topology instance")
                        .atom_ids()
                        .map(|atom| InstanceAtomId::new(id, atom))
                        .collect();
                    (Some(id), atoms)
                })
                .collect(),
            QeqGrouping::ConnectedComponents => {
                let mut selections = Vec::new();
                for (id, _) in topology.instances() {
                    for component in topology
                        .connected_components(id)
                        .expect("validated topology instance")
                    {
                        selections.push((Some(id), component));
                    }
                }
                selections
            }
        };
        let groups = selections
            .into_iter()
            .map(|(instance, atoms)| {
                let global_atoms = atoms
                    .iter()
                    .map(|atom| {
                        topology
                            .atom_index(*atom)
                            .expect("charge-group atom has a dense index")
                            .index()
                    })
                    .collect::<Vec<_>>();
                let formal_charge = atoms
                    .iter()
                    .map(|atom| {
                        f64::from(
                            topology
                                .atom(*atom)
                                .expect("charge-group atom is live")
                                .formal_charge,
                        )
                    })
                    .sum();
                Ok(PreparedChargeGroup {
                    instance,
                    system: system_from_selection(topology, reference, &atoms)?,
                    global_atoms,
                    formal_charge,
                })
            })
            .collect::<Result<Vec<_>, DreidingPrepareError>>()?;
        Ok(Self { whole, groups })
    }
}

#[derive(Debug)]
struct PreparedChargeGroup {
    instance: Option<MoleculeInstanceId>,
    system: System,
    global_atoms: Vec<usize>,
    formal_charge: f64,
}

fn validate_atoms(topology: &Topology) -> Result<(), DreidingPrepareError> {
    for (atom_id, atom) in topology.atoms() {
        let implicit = topology
            .implicit_hydrogens(atom_id)
            .expect("topology atom perception lookup");
        if atom.explicit_hydrogens != 0 || implicit.is_some_and(|count| count != 0) {
            return Err(DreidingPrepareError::CountedHydrogens {
                atom: atom_id,
                explicit: atom.explicit_hydrogens,
                implicit: implicit.unwrap_or(0),
            });
        }
        if implicit.is_none() && !atom.no_implicit_hydrogens {
            return Err(DreidingPrepareError::UnresolvedImplicitHydrogens { atom: atom_id });
        }
        if atom.radical.is_some() {
            return Err(DreidingPrepareError::RadicalAtom { atom: atom_id });
        }
    }
    Ok(())
}

fn validate_bonds(topology: &Topology) -> Result<(), DreidingPrepareError> {
    for (bond_id, bond) in topology.bonds() {
        let aromatic = topology
            .bond_is_aromatic(bond_id)
            .expect("topology bond perception lookup")
            .unwrap_or(false);
        forge_bond_order(bond_id, bond.order, aromatic)?;
    }
    Ok(())
}

fn system_from_selection(
    topology: &Topology,
    reference: ModelView<'_>,
    selected: &[InstanceAtomId],
) -> Result<System, DreidingPrepareError> {
    let local_by_global = selected
        .iter()
        .enumerate()
        .map(|(local, atom)| (*atom, local))
        .collect::<BTreeMap<_, _>>();
    let mut system = System::new();
    for &atom_id in selected {
        let atom = topology.atom(atom_id).expect("topology atom");
        let element = atom.element.symbol().parse().map_err(|_| {
            DreidingPrepareError::UnsupportedElement {
                atom: atom_id,
                symbol: atom.element.symbol().to_owned(),
            }
        })?;
        let point = reference
            .position(atom_id)
            .expect("complete reference positions")
            .into_value();
        system
            .atoms
            .push(ForgeAtom::new(element, [point.x, point.y, point.z]));
    }
    for (bond_id, bond) in topology.bonds() {
        let (a, b) = bond.endpoints();
        let a = InstanceAtomId::new(bond_id.molecule(), a);
        let b = InstanceAtomId::new(bond_id.molecule(), b);
        let (Some(&local_a), Some(&local_b)) = (local_by_global.get(&a), local_by_global.get(&b))
        else {
            continue;
        };
        system.bonds.push(ForgeBond::new(
            local_a,
            local_b,
            forge_bond_order(
                bond_id,
                bond.order,
                topology
                    .bond_is_aromatic(bond_id)
                    .expect("topology bond perception lookup")
                    .unwrap_or(false),
            )?,
        ));
    }
    Ok(system)
}

fn forge_bond_order(
    bond: InstanceBondId,
    order: BondOrder,
    aromatic: bool,
) -> Result<ForgeBondOrder, DreidingPrepareError> {
    match (aromatic, order) {
        (_, order @ (BondOrder::Zero | BondOrder::Quadruple | BondOrder::Dative)) => {
            Err(DreidingPrepareError::UnsupportedBondOrder { bond, order })
        }
        (false, BondOrder::Aromatic) | (true, BondOrder::Triple) => {
            Err(DreidingPrepareError::InconsistentAromaticBond { bond })
        }
        (true, BondOrder::Single | BondOrder::Double | BondOrder::Aromatic) => {
            Ok(ForgeBondOrder::Aromatic)
        }
        (false, BondOrder::Single) => Ok(ForgeBondOrder::Single),
        (false, BondOrder::Double) => Ok(ForgeBondOrder::Double),
        (false, BondOrder::Triple) => Ok(ForgeBondOrder::Triple),
    }
}

fn forge_system(
    system: &System,
    total_charge: f64,
    molecule: Option<MoleculeInstanceId>,
) -> Result<ForgedSystem, DreidingPrepareError> {
    let config = ForgeConfig {
        charge_method: ChargeMethod::Qeq(QeqConfig {
            total_charge,
            ..QeqConfig::default()
        }),
        bond_potential: BondPotentialType::Harmonic,
        angle_potential: AnglePotentialType::Cosine,
        vdw_potential: VdwPotentialType::LennardJones,
        ..ForgeConfig::default()
    };
    forge(system, &config).map_err(|error| DreidingPrepareError::Parameterization {
        molecule,
        message: error.to_string(),
    })
}

fn per_atom_types(forged: &ForgedSystem) -> Result<Vec<String>, DreidingPrepareError> {
    forged
        .atom_properties
        .iter()
        .map(|atom| {
            forged
                .atom_types
                .get(atom.type_idx)
                .cloned()
                .ok_or_else(|| DreidingPrepareError::InvalidPreparedData {
                    interaction: "atom type",
                    detail: format!("type index {} is unavailable", atom.type_idx),
                })
        })
        .collect()
}

fn prepare_bonds(forged: &ForgedSystem) -> Result<Vec<BondTerm>, DreidingPrepareError> {
    forged
        .potentials
        .bonds
        .iter()
        .map(|term| match *term {
            BondPotential::Harmonic { atoms, k_half, r0 } => {
                require_finite("bond", &[k_half, r0])?;
                Ok(BondTerm {
                    a: atoms.0,
                    b: atoms.1,
                    k_half: k_half * KCAL_TO_KJ,
                    r0,
                })
            }
            BondPotential::Morse { .. } => Err(DreidingPrepareError::InvalidPreparedData {
                interaction: "bond",
                detail: "forge returned Morse terms for the fixed harmonic configuration".into(),
            }),
        })
        .collect()
}

fn prepare_angles(forged: &ForgedSystem) -> Result<Vec<AngleTerm>, DreidingPrepareError> {
    forged
        .potentials
        .angles
        .iter()
        .map(|term| match *term {
            AnglePotential::CosineHarmonic {
                atoms,
                c_half,
                cos0,
            } => {
                require_finite("angle", &[c_half, cos0])?;
                Ok(AngleTerm::Harmonic {
                    atoms: [atoms.0, atoms.1, atoms.2],
                    c_half: c_half * KCAL_TO_KJ,
                    cos0,
                })
            }
            AnglePotential::CosineLinear { atoms, c } => {
                require_finite("angle", &[c])?;
                Ok(AngleTerm::Linear {
                    atoms: [atoms.0, atoms.1, atoms.2],
                    c: c * KCAL_TO_KJ,
                })
            }
            AnglePotential::ThetaHarmonic { .. } => {
                Err(DreidingPrepareError::InvalidPreparedData {
                    interaction: "angle",
                    detail:
                        "forge returned theta-harmonic terms for the fixed cosine configuration"
                            .into(),
                })
            }
        })
        .collect()
}

fn prepare_torsions(forged: &ForgedSystem) -> Result<Vec<TorsionTerm>, DreidingPrepareError> {
    forged
        .potentials
        .torsions
        .iter()
        .map(|term: &TorsionPotential| {
            require_finite("torsion", &[term.v_half, term.cos_n_phi0, term.sin_n_phi0])?;
            Ok(TorsionTerm {
                atoms: [term.atoms.0, term.atoms.1, term.atoms.2, term.atoms.3],
                v_half: term.v_half * KCAL_TO_KJ,
                n: term.n,
                cos_n_phi0: term.cos_n_phi0,
                sin_n_phi0: term.sin_n_phi0,
            })
        })
        .collect()
}

fn prepare_inversions(forged: &ForgedSystem) -> Result<Vec<InversionTerm>, DreidingPrepareError> {
    forged
        .potentials
        .inversions
        .iter()
        .map(|term| match *term {
            InversionPotential::Planar { atoms, c_half } => {
                require_finite("inversion", &[c_half])?;
                Ok(InversionTerm::Planar {
                    atoms: [atoms.0, atoms.1, atoms.2, atoms.3],
                    c_half: c_half * KCAL_TO_KJ,
                })
            }
            InversionPotential::Umbrella {
                atoms,
                c_half,
                cos_psi0,
            } => {
                require_finite("inversion", &[c_half, cos_psi0])?;
                Ok(InversionTerm::Umbrella {
                    atoms: [atoms.0, atoms.1, atoms.2, atoms.3],
                    c_half: c_half * KCAL_TO_KJ,
                    cos_psi0,
                })
            }
        })
        .collect()
}

fn prepare_nonbonded(
    topology: &Topology,
    forged: &ForgedSystem,
    charges: &[f64],
    exclusions: &BTreeSet<(usize, usize)>,
) -> Result<Vec<NonbondedTerm>, DreidingPrepareError> {
    let mut vdw = HashMap::new();
    for term in &forged.potentials.vdw_pairs {
        match *term {
            VdwPairPotential::LennardJones {
                type1_idx,
                type2_idx,
                d0,
                r0_sq,
            } => {
                require_finite("van der Waals", &[d0, r0_sq])?;
                vdw.insert(ordered_pair(type1_idx, type2_idx), (d0 * KCAL_TO_KJ, r0_sq));
            }
            VdwPairPotential::Buckingham { .. } => {
                return Err(DreidingPrepareError::InvalidPreparedData {
                    interaction: "van der Waals",
                    detail:
                        "forge returned Buckingham terms for the fixed Lennard-Jones configuration"
                            .into(),
                });
            }
        }
    }

    let mut terms = Vec::new();
    for first in 0..topology.atom_count() {
        for second in (first + 1)..topology.atom_count() {
            if exclusions.contains(&(first, second)) {
                continue;
            }
            let first_type = forged.atom_properties[first].type_idx;
            let second_type = forged.atom_properties[second].type_idx;
            let &(d0, r0_sq) = vdw.get(&ordered_pair(first_type, second_type)).ok_or(
                DreidingPrepareError::MissingVdwParameters {
                    first: topology.atom_ids()[first],
                    second: topology.atom_ids()[second],
                },
            )?;
            let coulomb = COULOMB_KJ_ANGSTROM_PER_MOL_E2 * charges[first] * charges[second];
            require_finite("electrostatic", &[coulomb])?;
            terms.push(NonbondedTerm {
                first,
                second,
                d0,
                r0_sq,
                coulomb,
            });
        }
    }
    Ok(terms)
}

fn prepare_hydrogen_bonds(
    forged: &ForgedSystem,
    adjacency: &[Vec<usize>],
    exclusions: &BTreeSet<(usize, usize)>,
) -> Result<Vec<HydrogenBondTerm>, DreidingPrepareError> {
    let mut parameters = HashMap::new();
    for term in &forged.potentials.h_bonds {
        let HBondPotential {
            donor_type_idx,
            hydrogen_type_idx,
            acceptor_type_idx,
            d_hb,
            r_hb_sq,
        } = *term;
        require_finite("hydrogen bond", &[d_hb, r_hb_sq])?;
        parameters.insert(
            (donor_type_idx, hydrogen_type_idx, acceptor_type_idx),
            (d_hb * KCAL_TO_KJ, r_hb_sq),
        );
    }

    let type_indices = forged
        .atom_properties
        .iter()
        .map(|atom| atom.type_idx)
        .collect::<Vec<_>>();
    let mut terms = Vec::new();
    let mut seen = BTreeSet::new();
    for hydrogen in 0..type_indices.len() {
        let Some(&donor) = adjacency[hydrogen].first() else {
            continue;
        };
        for acceptor in 0..type_indices.len() {
            if acceptor == donor || acceptor == hydrogen {
                continue;
            }
            if exclusions.contains(&ordered_pair(donor, acceptor)) {
                continue;
            }
            let key = (
                type_indices[donor],
                type_indices[hydrogen],
                type_indices[acceptor],
            );
            let Some(&(d_hb, r_hb_sq)) = parameters.get(&key) else {
                continue;
            };
            if seen.insert((donor, hydrogen, acceptor)) {
                terms.push(HydrogenBondTerm {
                    donor,
                    hydrogen,
                    acceptor,
                    d_hb,
                    r_hb_sq,
                });
            }
        }
    }
    Ok(terms)
}

fn adjacency(topology: &Topology) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); topology.atom_count()];
    for (bond_id, bond) in topology.bonds() {
        let (a, b) = bond.endpoints();
        let a = topology
            .atom_index(InstanceAtomId::new(bond_id.molecule(), a))
            .expect("topology bond endpoint must have a dense atom index")
            .index();
        let b = topology
            .atom_index(InstanceAtomId::new(bond_id.molecule(), b))
            .expect("topology bond endpoint must have a dense atom index")
            .index();
        adjacency[a].push(b);
        adjacency[b].push(a);
    }
    for neighbors in &mut adjacency {
        neighbors.sort_unstable();
    }
    adjacency
}

fn nonbonded_exclusions(adjacency: &[Vec<usize>]) -> BTreeSet<(usize, usize)> {
    let mut excluded = BTreeSet::new();
    for (center, neighbors) in adjacency.iter().enumerate() {
        for &neighbor in neighbors {
            excluded.insert(ordered_pair(center, neighbor));
        }
        for left in 0..neighbors.len() {
            for right in (left + 1)..neighbors.len() {
                excluded.insert(ordered_pair(neighbors[left], neighbors[right]));
            }
        }
    }
    excluded
}

fn ordered_pair(first: usize, second: usize) -> (usize, usize) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

fn require_finite(interaction: &'static str, values: &[f64]) -> Result<(), DreidingPrepareError> {
    if values.iter().all(|value| value.is_finite()) {
        Ok(())
    } else {
        Err(DreidingPrepareError::InvalidPreparedData {
            interaction,
            detail: "one or more values are non-finite".into(),
        })
    }
}

fn require_finite_slice(
    interaction: &'static str,
    values: &[f64],
) -> Result<(), DreidingPrepareError> {
    require_finite(interaction, values)
}

#[derive(Debug, Clone)]
pub(crate) struct BondTerm {
    pub(crate) a: usize,
    pub(crate) b: usize,
    pub(crate) k_half: f64,
    pub(crate) r0: f64,
}

#[derive(Debug, Clone)]
pub(crate) enum AngleTerm {
    Harmonic {
        atoms: [usize; 3],
        c_half: f64,
        cos0: f64,
    },
    Linear {
        atoms: [usize; 3],
        c: f64,
    },
}

#[derive(Debug, Clone)]
pub(crate) struct TorsionTerm {
    pub(crate) atoms: [usize; 4],
    pub(crate) v_half: f64,
    pub(crate) n: u8,
    pub(crate) cos_n_phi0: f64,
    pub(crate) sin_n_phi0: f64,
}

#[derive(Debug, Clone)]
pub(crate) enum InversionTerm {
    Planar {
        atoms: [usize; 4],
        c_half: f64,
    },
    Umbrella {
        atoms: [usize; 4],
        c_half: f64,
        cos_psi0: f64,
    },
}

#[derive(Debug, Clone)]
pub(crate) struct NonbondedTerm {
    pub(crate) first: usize,
    pub(crate) second: usize,
    pub(crate) d0: f64,
    pub(crate) r0_sq: f64,
    pub(crate) coulomb: f64,
}

#[derive(Debug, Clone)]
pub(crate) struct HydrogenBondTerm {
    pub(crate) donor: usize,
    pub(crate) hydrogen: usize,
    pub(crate) acceptor: usize,
    pub(crate) d_hb: f64,
    pub(crate) r_hb_sq: f64,
}
