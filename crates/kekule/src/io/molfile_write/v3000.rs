//! V3000 CTAB rendering of a prepared record.
use crate::algorithms::explicit_valence;
use crate::chemistry::SourceStereoBondMarkKind;
use crate::diagnostics::FormatError;
use crate::io::write_error::MOLFILE_WRITE;
use crate::molecule::{AtomRadical, BondOrder, ImplicitHydrogens};
use crate::stereo::StereoGroupKind;

use super::MolfileRecord;

/// CTfile lines hold at most 80 characters.
const MAX_LINE: usize = 80;
const PREFIX: &str = "M  V30 ";

pub(super) fn render(record: &MolfileRecord<'_>, title: &str) -> Result<String, FormatError> {
    let mut out = String::new();
    record.push_header(&mut out, title);
    out.push_str("  0  0  0  0  0  0            999 V3000\n");
    push_line(&mut out, "BEGIN CTAB");
    push_line(
        &mut out,
        &format!(
            "COUNTS {} {} 0 0 {}",
            record.atoms.len(),
            record.bonds.len(),
            u8::from(record.chiral)
        ),
    );
    push_line(&mut out, "BEGIN ATOM");
    for (index, record_atom) in (1u64..).zip(record.atoms.iter()) {
        let atom = record_atom.atom;
        let point = record_atom.position;
        let mut line = format!(
            "{index} {} {} {} {} {}",
            atom.element.symbol(),
            point.x,
            point.y,
            point.z,
            atom.atom_map.map_or(0, std::num::NonZero::get)
        );
        if atom.formal_charge != 0 {
            line.push_str(&format!(" CHG={}", atom.formal_charge));
        }
        if let Some(isotope) = atom.isotope {
            line.push_str(&format!(" MASS={isotope}"));
        }
        if let Some(radical) = atom.radical {
            line.push_str(&format!(" RAD={}", v3000_radical_code(radical)?));
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
                    line.push_str(" VAL=-1");
                } else {
                    line.push_str(&format!(" VAL={valence}"));
                }
            }
        }
        push_line(&mut out, &line);
    }
    push_line(&mut out, "END ATOM");
    push_line(&mut out, "BEGIN BOND");
    for (serial, bond) in (1u64..).zip(record.bonds.iter()) {
        let order_code = v3000_bond_code(bond.order)?;
        let mut line = format!("{serial} {order_code} {} {}", bond.from, bond.to);
        if let Some(cfg) = v3000_bond_cfg(bond.order, bond.stereo)? {
            line.push_str(&format!(" CFG={cfg}"));
        }
        push_line(&mut out, &line);
    }
    push_line(&mut out, "END BOND");
    if !record.stereo_groups.is_empty() {
        push_line(&mut out, "BEGIN COLLECTION");
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
                    return Err(
                        MOLFILE_WRITE.unsupported("V3000 cannot encode this stereo group kind")
                    )
                }
            };
            let mut line = format!("{kind} ATOMS=({}", atoms.len());
            for member in atoms {
                line.push_str(&format!(" {member}"));
            }
            line.push(')');
            push_line(&mut out, &line);
        }
        push_line(&mut out, "END COLLECTION");
    }
    push_line(&mut out, "END CTAB");
    out.push_str("M  END\n");
    Ok(out)
}

/// Writes one `M  V30` record, continuing it on following lines that end in
/// `-` when it exceeds the CTfile line length. Each break follows a space, so
/// no continuation starts with whitespace that readers trim; a single token
/// too long for one line is split inside.
fn push_line(out: &mut String, body: &str) {
    let room = MAX_LINE - PREFIX.len() - 1;
    let mut rest = body;
    while PREFIX.len() + rest.len() > MAX_LINE {
        // Readers strip one trailing `-` per line, so a token split inside
        // keeps any `-` of its own.
        let split = rest[..room].rfind(' ').map_or(room, |space| space + 1);
        out.push_str(PREFIX);
        out.push_str(&rest[..split]);
        out.push_str("-\n");
        rest = &rest[split..];
    }
    out.push_str(PREFIX);
    out.push_str(rest);
    out.push('\n');
}

fn v3000_bond_code(order: BondOrder) -> Result<u8, FormatError> {
    match order {
        BondOrder::Zero => Ok(0),
        BondOrder::Single => Ok(1),
        BondOrder::Double => Ok(2),
        BondOrder::Triple => Ok(3),
        BondOrder::Dative => Ok(9),
        BondOrder::Quadruple => {
            Err(MOLFILE_WRITE.unsupported("V3000 writer does not support quadruple bonds"))
        }
    }
}

fn v3000_bond_cfg(
    order: BondOrder,
    stereo: Option<SourceStereoBondMarkKind>,
) -> Result<Option<u8>, FormatError> {
    match (order, stereo) {
        (_, None) => Ok(None),
        (BondOrder::Double, Some(SourceStereoBondMarkKind::MolfileDoubleBondGeometry)) => Ok(None),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeUp)) => Ok(Some(1)),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeEither)) => Ok(Some(2)),
        (BondOrder::Single, Some(SourceStereoBondMarkKind::WedgeDown)) => Ok(Some(3)),
        (BondOrder::Double, Some(SourceStereoBondMarkKind::DoubleBondEither)) => Ok(Some(2)),
        _ => Err(MOLFILE_WRITE.unsupported("V3000 bond CFG is incompatible with the bond order")),
    }
}

fn v3000_radical_code(radical: AtomRadical) -> Result<u8, FormatError> {
    match (radical.electron_count(), radical.spin_multiplicity()) {
        (2, Some(1)) => Ok(1),
        (1, Some(2)) => Ok(2),
        (2, Some(3)) => Ok(3),
        _ => Err(MOLFILE_WRITE.unsupported(
            "V3000 writer requires an encodable radical electron count and explicit spin multiplicity",
        )),
    }
}
