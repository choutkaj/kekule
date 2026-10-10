//! Reads one small-molecule input into connected components keyed by source.

use kekule::core::Molecule;
use kekule::structure::{Model, Positions};
use kekule::topology::InstanceAtomId;
use kekule::units::{Quantity, ANGSTROM};
use kekule::{molfile, sdf, smiles};

use crate::Failure;

/// One connected molecule and the source position of each of its atoms.
pub struct Component {
    pub molecule: Molecule,
    /// Source atom index, by `AtomId::index()`.
    pub source: Vec<usize>,
    /// Coordinates in molecule atom order, when the source has them.
    pub positions: Option<Positions>,
}

impl Component {
    pub fn atom(&self, atom: kekule::core::AtomId) -> usize {
        self.source[atom.index()]
    }
}

/// Every component of one input record, in topology instance order.
pub fn read(format: &str, text: &str) -> Result<Vec<Component>, Failure> {
    match format {
        "smiles" => read_smiles(text),
        "mol" => {
            let document = molfile::parse_str(text).map_err(|e| Failure::new("parse", e))?;
            let interpretation = document.interpret().map_err(|e| Failure::new("parse", e))?;
            from_model(interpretation.model())
        }
        "sdf" => {
            let document = sdf::parse_str(text).map_err(|e| Failure::new("parse", e))?;
            let [record] = document.records() else {
                return Err(Failure::new(
                    "request",
                    format!(
                        "expected one SDF record, found {}",
                        document.records().len()
                    ),
                ));
            };
            let interpretation = record.interpret().map_err(|e| Failure::new("parse", e))?;
            from_model(interpretation.model())
        }
        other => Err(Failure::new("request", format!("unknown format {other:?}"))),
    }
}

fn read_smiles(text: &str) -> Result<Vec<Component>, Failure> {
    let line = text.lines().next().unwrap_or_default().trim();
    let document = smiles::parse_str(line).map_err(|e| Failure::new("parse", e))?;
    let interpretation = document.interpret().map_err(|e| Failure::new("parse", e))?;
    Ok(interpretation
        .components()
        .iter()
        .map(|component| {
            let molecule = component.molecule().clone();
            let mut source = vec![0; molecule.atom_count()];
            for mapping in component.report().atom_mappings() {
                source[mapping.atom().index()] = mapping.source_index();
            }
            Component {
                molecule,
                source,
                positions: None,
            }
        })
        .collect())
}

/// Molfile and SDF interpretation keeps source atom-block order as the dense
/// topology order, so a dense atom index is the source row.
fn from_model(model: &Model) -> Result<Vec<Component>, Failure> {
    let topology = model.topology();
    topology
        .molecules()
        .map(|instance| {
            let molecule = instance.molecule().clone();
            let mut source = Vec::with_capacity(molecule.atom_count());
            let mut points = Vec::with_capacity(molecule.atom_count());
            for atom in molecule.atom_ids() {
                let qualified = InstanceAtomId::new(instance.id(), atom);
                let dense = topology
                    .atom_index(qualified)
                    .ok_or_else(|| Failure::new("parse", "atom missing from dense layout"))?;
                source.push(dense.index());
                let position = model
                    .position(qualified)
                    .map_err(|e| Failure::new("parse", e))?;
                points.push(
                    position
                        .value_in(ANGSTROM)
                        .map_err(|e| Failure::new("parse", e))?,
                );
            }
            let positions = Positions::new(Quantity::new(points, ANGSTROM))
                .map_err(|e| Failure::new("parse", e))?;
            Ok(Component {
                molecule,
                source,
                positions: Some(positions),
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const INTERLEAVED: &str = "\n  test\n\n  3  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.0000    0.0000    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0\n    1.4000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  3  1  0\nM  CHG  1   2   1\nM  END\n";

    #[test]
    fn molfile_components_keep_their_atom_block_rows() {
        let components = read("mol", INTERLEAVED).unwrap();
        let rows = components
            .iter()
            .map(|c| {
                let symbols = c
                    .molecule
                    .atoms()
                    .map(|(id, atom)| (c.atom(id), atom.element.symbol()))
                    .collect::<Vec<_>>();
                (symbols, c.positions.as_ref().map(Positions::len))
            })
            .collect::<Vec<_>>();
        assert_eq!(
            rows,
            [
                (vec![(0, "C"), (2, "O")], Some(2)),
                (vec![(1, "Na")], Some(1)),
            ]
        );
    }
}
