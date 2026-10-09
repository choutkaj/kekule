//! Molfile and SDF output re-reads to the same chemistry in the same atom order.

use kekule::core::Molecule;
use kekule::geometry::Point3;
use kekule::molfile::{self, MolfileWriteOptions, MolfileWriteVersion};
use kekule::sdf;
use kekule::structure::{Model, Positions};
use kekule::units::{Quantity, ANGSTROM};

use crate::corpus::{assert_invariant, molecules};
use crate::support::cip_labels;

/// Element, isotope, charge, radical and total hydrogens of every atom, and
/// every bond as sorted endpoints with its localized order, in storage order.
fn chemistry(molecule: &Molecule) -> String {
    let atoms = molecule
        .atom_ids()
        .map(|id| {
            let atom = molecule.atom(id).unwrap();
            let hydrogens = molecule.total_hydrogens(id).unwrap();
            format!(
                "{}:{:?}:{}:{:?}:{hydrogens:?}",
                atom.element.symbol(),
                atom.isotope,
                atom.formal_charge,
                atom.radical
            )
        })
        .collect::<Vec<_>>();
    let mut bonds = molecule
        .bonds()
        .map(|(_, bond)| {
            let (a, b) = (bond.a().index(), bond.b().index());
            format!("{}-{}:{:?}", a.min(b), a.max(b), bond.order)
        })
        .collect::<Vec<_>>();
    bonds.sort();
    format!("{atoms:?} {bonds:?}")
}

/// Molecules read from Molfile or SDF keep their source coordinates; others
/// are written with zero coordinates unless they carry stereo that would need
/// a drawing. Two documented format limits are asserted rather than skipped:
/// Molfile cannot encode a radical without its spin multiplicity, and V2000
/// cannot encode coordinates wider than its fixed-width atom field.
#[test]
fn molfile_and_sdf_round_trips_preserve_chemistry_and_stereo() {
    let mut failures = Vec::new();
    let mut checked = 0;
    for sample in molecules() {
        let mut fail = |message: String| failures.push((sample.label.clone(), message));
        let model = match &sample.model {
            Some(model) => model.clone(),
            None if sample.molecule.stereo_elements().next().is_none() => {
                let zeros = vec![Point3::origin(); sample.molecule.atom_count()];
                let positions = Positions::new(Quantity::new(zeros, ANGSTROM)).unwrap();
                Model::from_molecule(sample.molecule.clone(), &positions).unwrap()
            }
            None => continue,
        };
        let unspecified_spin = sample.molecule.atoms().any(|(_, atom)| {
            atom.radical
                .is_some_and(|r| r.spin_multiplicity().is_none())
        });
        // Molfiles store angstroms in `%10.4f` columns.
        let angstroms = model
            .positions()
            .values()
            .map(|points| points.to_vec())
            .value_in(ANGSTROM)
            .unwrap();
        let v2000_fits = angstroms.iter().all(|point| {
            [point.x, point.y, point.z]
                .into_iter()
                .all(|coordinate| format!("{coordinate:.4}").len() <= 10)
        });
        let mut outputs = Vec::new();
        for (format, written) in [
            (
                "V2000",
                molfile::write(
                    &model,
                    MolfileWriteOptions {
                        version: MolfileWriteVersion::V2000,
                    },
                )
                .map_err(|error| error.to_string()),
            ),
            (
                "V3000",
                molfile::write(
                    &model,
                    MolfileWriteOptions {
                        version: MolfileWriteVersion::V3000,
                    },
                )
                .map_err(|error| error.to_string()),
            ),
            (
                "SDF",
                sdf::write([&model], Default::default()).map_err(|error| error.to_string()),
            ),
        ] {
            let must_reject = unspecified_spin || (format == "V2000" && !v2000_fits);
            match (written, must_reject) {
                (Ok(text), false) => outputs.push((format, text)),
                (Err(_), true) => {}
                (Ok(_), true) => fail(format!("{format}: writer accepted unencodable content")),
                (Err(error), false) => fail(format!("{format}: {error}")),
            }
        }
        if outputs.is_empty() {
            continue;
        }
        checked += 1;
        let mut expected = sample.perceived();
        let expected_chemistry = chemistry(&expected);
        let expected_labels = cip_labels(&mut expected);
        for (format, text) in outputs {
            let reread = if format == "SDF" {
                sdf::parse_str(&text)
                    .map_err(|e| e.to_string())
                    .and_then(|document| {
                        document.records()[0].to_model().map_err(|e| e.to_string())
                    })
            } else {
                molfile::parse_str(&text)
                    .map_err(|e| e.to_string())
                    .and_then(|document| document.to_model().map_err(|e| e.to_string()))
            };
            let reread = match reread {
                Ok(model) => model,
                Err(error) => {
                    fail(format!("{format}: re-read failed: {error}"));
                    continue;
                }
            };
            let mut molecules = reread
                .topology()
                .molecules()
                .map(|instance| instance.molecule().clone())
                .collect::<Vec<_>>();
            if molecules.len() != 1 {
                fail(format!("{format}: {} components", molecules.len()));
                continue;
            }
            let mut actual = molecules.pop().unwrap();
            actual.perceive().unwrap();
            if chemistry(&actual) != expected_chemistry {
                fail(format!(
                    "{format}: chemistry changed
  {}
  {expected_chemistry}",
                    chemistry(&actual)
                ));
            }
            let labels = cip_labels(&mut actual);
            if labels != expected_labels {
                fail(format!("{format}: CIP {labels:?} != {expected_labels:?}"));
            }
        }
    }
    assert!(checked > 100, "only {checked} molecules were written");
    assert_invariant("molfile_round_trip", failures);
}
