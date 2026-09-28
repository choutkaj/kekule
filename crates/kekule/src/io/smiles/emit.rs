//! Structured SMILES emission; CX atom indices are assigned after traversal.
use std::collections::{BTreeMap, BTreeSet};

use crate::core::{AtomId, Molecule, StereoElementKind, StereoGroupKind};
use crate::io::MolWriteError;

#[derive(Debug, Clone, Default)]
pub(super) struct Emission {
    pub text: String,
    pub atoms: usize,
    groups: Vec<(StereoGroupKind, Vec<usize>)>,
    ungrouped: Vec<usize>,
    radicals: Vec<(u8, usize)>,
}

impl Emission {
    pub fn new(mol: &Molecule, text: String, order: &[AtomId]) -> Self {
        let indices: BTreeMap<_, _> = order.iter().enumerate().map(|(i, &a)| (a, i)).collect();
        let centers: BTreeMap<_, _> = mol
            .stereo_elements()
            .filter_map(|(id, e)| {
                if let StereoElementKind::Tetrahedral(s) = &e.kind {
                    indices.get(&s.center).map(|&i| (id, i))
                } else {
                    None
                }
            })
            .collect();
        let groups = mol
            .stereo_groups()
            .filter_map(|(_, g)| {
                let mut members: Vec<_> = g
                    .members
                    .iter()
                    .filter_map(|id| centers.get(id).copied())
                    .collect();
                members.sort_unstable();
                (!members.is_empty()).then_some((g.kind, members))
            })
            .collect();
        let ungrouped = mol
            .stereo_elements()
            .filter(|(_, e)| e.group.is_none())
            .filter_map(|(id, _)| centers.get(&id).copied())
            .collect();
        Self {
            text,
            atoms: order.len(),
            groups,
            ungrouped,
            radicals: order
                .iter()
                .enumerate()
                .filter_map(|(i, atom)| {
                    let radical = mol
                        .atom(*atom)
                        .expect("emission visits a live atom")
                        .radical?;
                    // Explicit spin remains unsupported by the writer. Higher
                    // occupancies have no CX code and remain encoded by brackets.
                    let code = match radical.electron_count() {
                        1 => 1,
                        2 => 2,
                        3 => 5,
                        _ => return None,
                    };
                    Some((code, i))
                })
                .collect(),
        }
    }

    pub fn join(parts: Vec<Self>) -> Self {
        let mut result = Self::default();
        for mut part in parts {
            if !result.text.is_empty() {
                result.text.push('.');
            }
            result.text.push_str(&part.text);
            for (_, members) in &mut part.groups {
                for member in members {
                    *member += result.atoms;
                }
            }
            result
                .ungrouped
                .extend(part.ungrouped.into_iter().map(|i| i + result.atoms));
            result.groups.extend(part.groups);
            result.radicals.extend(
                part.radicals
                    .into_iter()
                    .map(|(code, i)| (code, i + result.atoms)),
            );
            result.atoms += part.atoms;
        }
        result
    }

    pub fn render(&self) -> Result<String, MolWriteError> {
        let mut absolute = BTreeSet::new();
        let mut and = Vec::new();
        let mut or = Vec::new();
        let mut relative = 0;
        for (kind, members) in &self.groups {
            match kind {
                StereoGroupKind::Absolute => absolute.extend(members.iter().copied()),
                StereoGroupKind::And => and.push(members.clone()),
                StereoGroupKind::Or => or.push(members.clone()),
                StereoGroupKind::Relative => relative += 1,
                StereoGroupKind::Racemic => return Err(MolWriteError::new(
                    "CXSMILES cannot encode a quantitative racemic group without losing its distinction from AND")),
            }
        }
        if relative > 1 {
            return Err(MolWriteError::new(
                "CXSMILES r cannot encode multiple independent relative stereo groups",
            ));
        }
        // The record-level r flag applies to every otherwise ungrouped center.
        // Shield known absolute centers so it cannot weaken their configuration.
        if relative != 0 {
            absolute.extend(self.ungrouped.iter().copied());
        }
        let list = |members: &[usize]| {
            members
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",")
        };
        let mut fields = Vec::new();
        let mut radicals = BTreeMap::<u8, Vec<usize>>::new();
        for &(code, atom) in &self.radicals {
            radicals.entry(code).or_default().push(atom);
        }
        for (code, atoms) in radicals {
            fields.push(format!("^{code}:{}", list(&atoms)));
        }
        if !absolute.is_empty() {
            fields.push(format!(
                "a:{}",
                list(&absolute.into_iter().collect::<Vec<_>>())
            ));
        }
        and.sort();
        or.sort();
        for (tag, groups) in [("o", or), ("&", and)] {
            for (number, members) in groups.iter().enumerate() {
                fields.push(format!("{tag}{}:{}", number + 1, list(members)));
            }
        }
        if relative != 0 {
            fields.push("r".to_owned());
        }
        if fields.is_empty() {
            Ok(self.text.clone())
        } else {
            Ok(format!("{} |{}|", self.text, fields.join(",")))
        }
    }
}

pub(super) fn validate_groups(mol: &Molecule) -> Result<(), MolWriteError> {
    let mut relative = 0;
    for (_, group) in mol.stereo_groups() {
        if group.kind == StereoGroupKind::Racemic {
            return Err(MolWriteError::new("CXSMILES cannot encode a quantitative racemic group without losing its distinction from AND"));
        }
        relative += usize::from(group.kind == StereoGroupKind::Relative);
        for member in &group.members {
            if !matches!(
                mol.stereo_element(*member).map(|s| &s.kind),
                Ok(StereoElementKind::Tetrahedral(_))
            ) {
                return Err(MolWriteError::new("CXSMILES group export requires tetrahedral members; this stereo family has no supported SMILES encoding"));
            }
        }
    }
    if relative > 1 {
        return Err(MolWriteError::new(
            "CXSMILES r cannot encode multiple independent relative stereo groups",
        ));
    }
    Ok(())
}
