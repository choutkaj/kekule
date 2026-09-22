use super::atom_hydrogen_count;
use crate::core::*;

/// Local tetrahedral carrier geometry. Bond multiplicity changes electron
/// counting, not the number of ligand directions (for example a phosphoryl O).
pub(super) fn tetrahedral_carriers(mol: &Molecule, center: AtomId) -> Option<Vec<StereoCarrier>> {
    let atom = mol.atom(center).ok()?;
    if atom.element.symbol() == "H" {
        return None;
    }
    let mut carriers = Vec::new();
    let mut valence = usize::from(atom_hydrogen_count(mol, center));
    for (_, bond) in mol.incident_bonds(center).ok()? {
        valence += match bond.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            // These coordination environments are outside this tetrahedral model.
            _ => return None,
        };
        carriers.push(StereoCarrier::Atom(bond.other_atom(center)));
    }
    let hydrogens = atom_hydrogen_count(mol, center);
    if hydrogens > 1 {
        return None;
    }
    carriers.sort_by_key(|carrier| carrier.canonical_order_key());
    if hydrogens == 1 {
        carriers.push(StereoCarrier::ImplicitHydrogen);
    }
    if carriers.len() == 3 {
        let lone_pair = match atom.element.symbol() {
            // Trigonal-pyramidal phosphines/arsines may have a hydrogen ligand.
            "P" | "As" => true,
            "N" => super::nitrogen::lone_pair_is_stereogenic(mol, center),
            "S" | "Se" => {
                hydrogens == 0 && (valence == 4 || (valence == 3 && atom.formal_charge == 1))
            }
            _ => false,
        };
        if lone_pair {
            carriers.push(StereoCarrier::ImplicitLonePair);
        }
    }
    (carriers.len() == 4).then_some(carriers)
}

pub(super) fn unclassified_tetrahedral_geometry(mol: &Molecule, center: AtomId) -> bool {
    tetrahedral_carriers(mol, center).is_none()
        && atom_hydrogen_count(mol, center) <= 1
        && !super::nitrogen::classified(mol, center)
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct AtomIdentity {
    element: u8,
    isotope: Option<u16>,
    charge: i8,
    radical: Option<(u8, Option<u8>)>,
}

impl AtomIdentity {
    pub(super) fn of(atom: &Atom) -> Self {
        Self {
            element: atom.element.atomic_number(),
            isotope: atom.isotope,
            charge: atom.formal_charge,
            radical: atom
                .radical
                .map(|radical| (radical.electron_count(), radical.spin_multiplicity())),
        }
    }
    pub(super) fn hydrogen() -> Self {
        Self {
            element: 1,
            isotope: None,
            charge: 0,
            radical: None,
        }
    }
}

#[derive(PartialEq, Eq)]
struct TerminalLigand {
    atom: AtomIdentity,
    bond: BondOrder,
    hydrogens: Vec<AtomIdentity>,
}

/// Proves equivalence only for terminal atoms/groups whose other neighbors are
/// terminal hydrogens. General branch and ring equivalence needs a separate
/// symmetry analysis. Atom maps are annotations, not chemical distinctions.
pub(super) fn has_repeated_terminal_ligands(
    mol: &Molecule,
    center: AtomId,
    carriers: &[StereoCarrier],
) -> bool {
    let mut known = Vec::new();
    for carrier in carriers {
        let identity = match carrier {
            StereoCarrier::ImplicitHydrogen => Some(TerminalLigand {
                atom: AtomIdentity::hydrogen(),
                bond: BondOrder::Single,
                hydrogens: Vec::new(),
            }),
            StereoCarrier::Atom(atom) => terminal_ligand(mol, center, *atom),
            StereoCarrier::ImplicitLonePair => None,
        };
        if let Some(identity) = identity {
            if known.contains(&identity) {
                return true;
            }
            known.push(identity);
        }
    }
    false
}

fn terminal_ligand(mol: &Molecule, center: AtomId, carrier: AtomId) -> Option<TerminalLigand> {
    let atom = mol.atom(carrier).ok()?;
    let count = mol.implicit_hydrogens(carrier).ok()??;
    let mut hydrogens: Vec<_> = (0..count).map(|_| AtomIdentity::hydrogen()).collect();
    let mut attachment = None;
    for (_, bond) in mol.incident_bonds(carrier).ok()? {
        let neighbor = bond.other_atom(carrier);
        if neighbor == center {
            attachment = Some(bond.order);
            continue;
        }
        let hydrogen = mol.atom(neighbor).ok()?;
        if hydrogen.element.symbol() != "H"
            || bond.order != BondOrder::Single
            || mol.neighbors(neighbor).ok()?.count() != 1
            || atom_hydrogen_count(mol, neighbor) != 0
        {
            return None;
        }
        hydrogens.push(AtomIdentity::of(hydrogen));
    }
    hydrogens.sort_unstable();
    // Isotopically distinct hydrogen ligands can make this terminal group
    // stereogenic itself. Its relationship to another such group then belongs
    // to the stereo dependency analysis, not this constitutional proof.
    if hydrogens.len() >= 2 && hydrogens.windows(2).all(|pair| pair[0] != pair[1]) {
        return None;
    }
    Some(TerminalLigand {
        atom: AtomIdentity::of(atom),
        bond: attachment?,
        hydrogens,
    })
}
