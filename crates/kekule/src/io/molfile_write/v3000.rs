//! V3000 CTAB rendering of a prepared record.
use crate::algorithms::explicit_valence;
use crate::chemistry::SourceStereoBondMarkKind;
use crate::core::{AtomRadical, BondOrder, ImplicitHydrogens, StereoGroupKind};
use crate::io::MoleculeWriteError;

use super::MolfileRecord;

pub(super) fn render(
    record: &MolfileRecord<'_>,
    title: &str,
) -> Result<String, MoleculeWriteError> {
    let mut out = String::new();
    record.push_header(&mut out, title);
    out.push_str("  0  0  0  0  0  0            999 V3000\n");
    out.push_str("M  V30 BEGIN CTAB\n");
    out.push_str(&format!(
        "M  V30 COUNTS {} {} 0 0 0\n",
        record.atoms.len(),
        record.bonds.len()
    ));
    out.push_str("M  V30 BEGIN ATOM\n");
    for (index, record_atom) in (1u64..).zip(record.atoms.iter()) {
        let atom = record_atom.atom;
        let point = record_atom.position;
        out.push_str(&format!(
            "M  V30 {index} {} {} {} {} {}",
            atom.element.symbol(),
            point.x,
            point.y,
            point.z,
            atom.atom_map.unwrap_or(0)
        ));
        if atom.formal_charge != 0 {
            out.push_str(&format!(" CHG={}", atom.formal_charge));
        }
        if let Some(isotope) = atom.isotope {
            out.push_str(&format!(" MASS={isotope}"));
        }
        if let Some(radical) = atom.radical {
            out.push_str(&format!(" RAD={}", v3000_radical_code(radical)?));
        }
        match atom.hydrogens {
            ImplicitHydrogens::Inferred => {}
            ImplicitHydrogens::Fixed(explicit) => {
                // HCOUNT is a query constraint in CTfile, not a molecular
                // hydrogen declaration. VAL preserves the fixed total valence
                // without turning the emitted atom into a query atom.
                let valence =
                    explicit_valence(record_atom.molecule, record_atom.id) + usize::from(explicit);
                // VAL=0 means unspecified, whereas -1 explicitly means zero.
                if valence == 0 {
                    out.push_str(" VAL=-1");
                } else {
                    out.push_str(&format!(" VAL={valence}"));
                }
            }
        }
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");
    out.push_str("M  V30 BEGIN BOND\n");
    for (serial, bond) in (1u64..).zip(record.bonds.iter()) {
        let order_code = v3000_bond_code(bond.order)?;
        out.push_str(&format!(
            "M  V30 {serial} {order_code} {} {}",
            bond.from, bond.to
        ));
        if let Some(cfg) = v3000_bond_cfg(bond.order, bond.stereo)? {
            out.push_str(&format!(" CFG={cfg}"));
        }
        out.push('\n');
    }
    out.push_str("M  V30 END BOND\n");
    if !record.stereo_groups.is_empty() {
        out.push_str("M  V30 BEGIN COLLECTION\n");
        let mut absolute_written = false;
        for (index, group) in record.stereo_groups.iter().enumerate() {
            if group.kind == StereoGroupKind::Absolute && absolute_written {
                continue;
            }
            let atoms = if group.kind == StereoGroupKind::Absolute {
                absolute_written = true;
                record
                    .stereo_groups
                    .iter()
                    .filter(|group| group.kind == StereoGroupKind::Absolute)
                    .flat_map(|group| group.atoms.iter().copied())
                    .collect::<Vec<_>>()
            } else {
                group.atoms.clone()
            };
            let kind = match group.kind {
                StereoGroupKind::Absolute => "MDLV30/STEABS".to_owned(),
                StereoGroupKind::Or => format!("MDLV30/STEREL{}", index + 1),
                StereoGroupKind::And => format!("MDLV30/STERAC{}", index + 1),
                _ => {
                    return Err(MoleculeWriteError::unsupported(
                        "V3000 cannot encode this stereo group kind",
                    ))
                }
            };
            out.push_str(&format!("M  V30 {kind}"));
            out.push_str(&format!(" ATOMS=({}", atoms.len()));
            for member in atoms {
                out.push_str(&format!(" {member}"));
            }
            out.push(')');
            out.push('\n');
        }
        out.push_str("M  V30 END COLLECTION\n");
    }
    out.push_str("M  V30 END CTAB\n");
    out.push_str("M  END\n");
    Ok(out)
}

fn v3000_bond_code(order: BondOrder) -> Result<u8, MoleculeWriteError> {
    match order {
        BondOrder::Zero => Ok(0),
        BondOrder::Single => Ok(1),
        BondOrder::Double => Ok(2),
        BondOrder::Triple => Ok(3),
        BondOrder::Dative => Ok(9),
        BondOrder::Quadruple => Err(MoleculeWriteError::unsupported(
            "V3000 writer does not support quadruple bonds",
        )),
    }
}

fn v3000_bond_cfg(
    order: BondOrder,
    stereo: Option<SourceStereoBondMarkKind>,
) -> Result<Option<u8>, MoleculeWriteError> {
    match (order, stereo) {
        (_, None) => Ok(None),
        (BondOrder::Double, Some(SourceStereoBondMarkKind::MolfileDoubleBondGeometry)) => Ok(None),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeUp)) => Ok(Some(1)),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeEither)) => Ok(Some(2)),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeDown)) => Ok(Some(3)),
        (BondOrder::Double, Some(SourceStereoBondMarkKind::DoubleBondEither)) => Ok(Some(2)),
        _ => Err(MoleculeWriteError::unsupported(
            "V3000 bond CFG is incompatible with the bond order",
        )),
    }
}

fn v3000_radical_code(radical: AtomRadical) -> Result<u8, MoleculeWriteError> {
    match (radical.electron_count(), radical.spin_multiplicity()) {
        (2, Some(1)) => Ok(1),
        (1, Some(2)) => Ok(2),
        (2, Some(3)) => Ok(3),
        _ => Err(MoleculeWriteError::unsupported(
            "V3000 writer requires an encodable radical electron count and explicit spin multiplicity",
        )),
    }
}
