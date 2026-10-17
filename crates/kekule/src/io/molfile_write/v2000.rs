//! V2000 CTAB rendering of a prepared record.
use crate::algorithms::explicit_valence;
use crate::chemistry::SourceStereoBondMarkKind;
use crate::diagnostics::FormatError;
use crate::io::v2000::{V2000_MAX_ATOMS, V2000_MAX_BONDS};
use crate::io::write_error::MOLFILE_WRITE;
use crate::molecule::{Atom, AtomId, AtomRadical, BondOrder, ImplicitHydrogens, Molecule};

use super::MolfileRecord;

/// Coordinates occupy fixed 10-column fields with four decimal places.
const COORDINATE_DECIMALS: usize = 4;

/// The value a coordinate takes once written, so stereo projection and the
/// header's dimensional code can judge the geometry readers will see.
pub(super) fn emitted_coordinate(value: f64) -> f64 {
    format!("{value:.COORDINATE_DECIMALS$}")
        .parse()
        .expect("a formatted coordinate is a valid floating-point number")
}

/// Rejects records whose size or enhanced stereo V2000 cannot hold. It needs
/// only counts, so a caller can test a source before preparing a record.
pub(super) fn check_capacity(
    atoms: usize,
    bonds: usize,
    has_stereo_groups: bool,
) -> Result<(), FormatError> {
    if has_stereo_groups {
        return Err(
            MOLFILE_WRITE.unsupported("V2000 cannot encode enhanced stereo groups; use V3000")
        );
    }
    if atoms > V2000_MAX_ATOMS || bonds > V2000_MAX_BONDS {
        return Err(
            MOLFILE_WRITE.unsupported("V2000 writer supports at most 999 atoms and 999 bonds")
        );
    }
    Ok(())
}

pub(super) fn render(record: &MolfileRecord<'_>, title: &str) -> Result<String, FormatError> {
    check_capacity(
        record.atoms.len(),
        record.bonds.len(),
        !record.stereo_groups.is_empty(),
    )?;

    let mut out = String::new();
    record.push_header(&mut out, title);
    out.push_str(&format!(
        "{:>3}{:>3}  0  0{:>3}  0            999 V2000\n",
        record.atoms.len(),
        record.bonds.len(),
        u8::from(record.chiral)
    ));

    for record_atom in &record.atoms {
        let atom = record_atom.atom;
        let point = record_atom.position;
        if atom.atom_map.is_some_and(|atom_map| atom_map > 999) {
            return Err(MOLFILE_WRITE.unsupported(format!(
                "V2000 cannot encode atom-map number {} above 999",
                atom.atom_map.expect("atom map was checked")
            )));
        }
        if [point.x, point.y, point.z]
            .into_iter()
            .any(|coordinate| format!("{coordinate:.COORDINATE_DECIMALS$}").len() > 10)
        {
            return Err(MOLFILE_WRITE
                .unsupported("V2000 cannot encode a coordinate in its fixed-width atom field"));
        }
        let valence_code =
            v2000_valence_code(record_atom.molecule, record_atom.id, record_atom.atom)?;
        out.push_str(&format!(
            "{:>10.COORDINATE_DECIMALS$}{:>10.COORDINATE_DECIMALS$}{:>10.COORDINATE_DECIMALS$} {:<3}{:>2}{:>3}  0  0  0{:>3}  0  0  0{:>3}  0  0\n",
            point.x,
            point.y,
            point.z,
            atom.element.symbol(),
            0,
            v2000_charge_code(atom.formal_charge),
            valence_code,
            atom.atom_map.unwrap_or(0)
        ));
    }

    for bond in &record.bonds {
        let order_code = v2000_bond_code(bond.order)?;
        let stereo_code = v2000_bond_stereo_code(bond.order, bond.stereo)?;
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}{:>3}  0  0  0\n",
            bond.from, bond.to, order_code, stereo_code
        ));
    }

    push_m_record(
        &mut out,
        "CHG",
        record
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(index, record_atom)| {
                (record_atom.atom.formal_charge != 0).then_some((
                    (index + 1) as u64,
                    i32::from(record_atom.atom.formal_charge),
                ))
            })
            .collect(),
    );
    push_m_record(
        &mut out,
        "ISO",
        record
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(index, record_atom)| {
                record_atom
                    .atom
                    .isotope
                    .map(|isotope| ((index + 1) as u64, i32::from(isotope)))
            })
            .collect(),
    );
    let mut radical_records = Vec::new();
    for (index, record_atom) in record.atoms.iter().enumerate() {
        let Some(radical) = record_atom.atom.radical else {
            continue;
        };
        radical_records.push(((index + 1) as u64, v2000_radical_code(radical)?));
    }
    push_m_record(&mut out, "RAD", radical_records);
    out.push_str("M  END\n");
    Ok(out)
}

fn v2000_charge_code(charge: i8) -> i8 {
    match charge {
        3 => 1,
        2 => 2,
        1 => 3,
        -1 => 5,
        -2 => 6,
        -3 => 7,
        _ => 0,
    }
}

fn v2000_valence_code(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> Result<u8, FormatError> {
    let fixed_hydrogens = match atom.hydrogens {
        ImplicitHydrogens::Inferred => return Ok(0),
        ImplicitHydrogens::Fixed(explicit) => explicit,
    };
    let valence = explicit_valence(mol, atom_id) + usize::from(fixed_hydrogens);
    match valence {
        0 => Ok(15),
        1..=14 => Ok(u8::try_from(valence).expect("range checked")),
        _ => Err(MOLFILE_WRITE.unsupported(format!(
            "V2000 cannot encode explicit valence {valence} for atom {}",
            atom_id.index()
        ))),
    }
}

fn v2000_bond_code(order: BondOrder) -> Result<u8, FormatError> {
    match order {
        BondOrder::Zero => Ok(0),
        BondOrder::Single => Ok(1),
        BondOrder::Double => Ok(2),
        BondOrder::Triple => Ok(3),
        BondOrder::Dative => Ok(9),
        BondOrder::Quadruple => {
            Err(MOLFILE_WRITE.unsupported("V2000 writer does not support quadruple bonds"))
        }
    }
}

fn v2000_bond_stereo_code(
    order: BondOrder,
    stereo: Option<SourceStereoBondMarkKind>,
) -> Result<u8, FormatError> {
    match (order, stereo) {
        (BondOrder::Double, Some(SourceStereoBondMarkKind::MolfileDoubleBondGeometry)) => Ok(0),
        (_, None) => Ok(0),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeUp)) => Ok(1),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeEither)) => Ok(4),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeDown)) => Ok(6),
        (BondOrder::Double, Some(SourceStereoBondMarkKind::DoubleBondEither)) => Ok(3),
        _ => {
            Err(MOLFILE_WRITE.unsupported("V2000 bond stereo is incompatible with the bond order"))
        }
    }
}

fn v2000_radical_code(radical: AtomRadical) -> Result<i32, FormatError> {
    match (radical.electron_count(), radical.spin_multiplicity()) {
        (2, Some(1)) => Ok(1),
        (1, Some(2)) => Ok(2),
        (2, Some(3)) => Ok(3),
        _ => Err(MOLFILE_WRITE.unsupported(
            "V2000 writer requires an encodable radical electron count and explicit spin multiplicity",
        )),
    }
}

fn push_m_record(out: &mut String, code: &str, pairs: Vec<(u64, i32)>) {
    for chunk in pairs.chunks(8) {
        if chunk.is_empty() {
            continue;
        }
        out.push_str(&format!("M  {code} {:>2}", chunk.len()));
        for (atom, value) in chunk {
            out.push_str(&format!("{atom:>4}{value:>4}"));
        }
        out.push('\n');
    }
}
