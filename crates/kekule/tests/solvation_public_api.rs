use kekule::{
    core::{Element, VanDerWaalsRadiusSource},
    geometry::{PeriodicCell, PeriodicGeometry, PeriodicGeometryError, Point3, Vector3},
    properties::{PropertyKey, PropertyValue},
    structure::{
        BoxShape, Model, NegativeIon, PeriodicBoxError, PeriodicBoxOptions, Positions, PositiveIon,
        SolvationChargeBasis, SolvationError, SolvationRadii, SolventOptions,
    },
    topology::{AtomSiteMetadata, MoleculeClass, ResidueClass},
    units::{Quantity, ANGSTROM, METER, MOLAR, MOLE, NANOMETER},
};
use std::sync::Arc;

fn model(smiles: &str) -> Model {
    let m = kekule::smiles::to_molecules(smiles).unwrap().remove(0);
    Model::from_molecule(&m, &Positions::zeros(m.atom_count())).unwrap()
}
fn boxed(smiles: &str) -> Model {
    let mut m = model(smiles);
    m.add_periodic_box(&PeriodicBoxOptions::Dimensions(Quantity::new(
        Vector3::new(3.0, 3.0, 3.0),
        NANOMETER,
    )))
    .unwrap();
    m
}
fn options() -> SolventOptions {
    SolventOptions {
        neutralize: false,
        ..Default::default()
    }
}
fn points(m: &Model) -> Vec<Point3> {
    m.positions().values().value().to_vec()
}
fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

#[test]
fn box_shapes_units_padding_and_snapshot_identity() {
    let mut m = model("CC");
    let atoms = m.topology().atom_ids().to_vec();
    m.set_atom_positions([
        (
            atoms[0],
            Quantity::new(Point3::new(20.0, 30.0, 40.0), NANOMETER),
        ),
        (
            atoms[1],
            Quantity::new(Point3::new(22.0, 30.0, 40.0), NANOMETER),
        ),
    ])
    .unwrap();
    let source = m.clone();
    // D=2, p=1, so width=3 (not D+2p).
    for (shape, ratio) in [
        (BoxShape::Cube, 1.0),
        (BoxShape::Dodecahedron, 2.0_f64.sqrt() / 2.0),
        (BoxShape::Octahedron, 4.0 / (3.0 * 3.0_f64.sqrt())),
    ] {
        m.add_periodic_box(&PeriodicBoxOptions::Padding {
            padding: Quantity::new(10.0, ANGSTROM),
            shape,
        })
        .unwrap();
        assert!((m.cell().unwrap().signed_volume() - 27.0 * ratio).abs() < 1e-12);
        assert_eq!(m.positions(), source.positions());
        assert!(Arc::ptr_eq(&m.shared_topology(), &source.shared_topology()));
        let geometry = PeriodicGeometry::new(*m.cell().unwrap()).unwrap();
        assert!((geometry.shortest_translation().unwrap().norm() - 3.0).abs() < 1e-12);
    }
    let mut single = model("[Na+]");
    single
        .add_periodic_box(&PeriodicBoxOptions::Padding {
            padding: Quantity::new(1.0, NANOMETER),
            shape: BoxShape::Cube,
        })
        .unwrap();
    assert_eq!(single.cell().unwrap().signed_volume(), 8.0);
    let cell = PeriodicCell::new(
        Quantity::new(
            [
                Vector3::new(0.0, 2.0, 0.0),
                Vector3::new(3.0, 0.0, 0.0),
                Vector3::new(0.1, 0.2, 4.0),
            ],
            NANOMETER,
        ),
        [true; 3],
    )
    .unwrap();
    m.add_periodic_box(&PeriodicBoxOptions::Cell(cell)).unwrap();
    assert_eq!(m.cell(), Some(&cell));
    m.add_periodic_box(&PeriodicBoxOptions::Dimensions(Quantity::new(
        Vector3::new(20.0, 30.0, 40.0),
        ANGSTROM,
    )))
    .unwrap();
    assert!((m.cell().unwrap().signed_volume() - 24.0).abs() < 1e-12);
}
#[test]
fn invalid_box_inputs_are_atomic() {
    let mut m = boxed("C");
    let before = m.clone();
    for padding in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(m
            .add_periodic_box(&PeriodicBoxOptions::Padding {
                padding: Quantity::new(padding, NANOMETER),
                shape: BoxShape::Cube
            })
            .is_err());
        assert_eq!(m, before);
    }
    assert!(m
        .add_periodic_box(&PeriodicBoxOptions::Padding {
            padding: Quantity::new(1.0, MOLAR),
            shape: BoxShape::Cube
        })
        .is_err());
    assert!(m
        .add_periodic_box(&PeriodicBoxOptions::Dimensions(Quantity::new(
            Vector3::new(0.0, 1.0, 1.0),
            NANOMETER
        )))
        .is_err());
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
        [true, false, true],
    )
    .unwrap();
    assert!(matches!(
        m.add_periodic_box(&PeriodicBoxOptions::Cell(cell)),
        Err(PeriodicBoxError::NotFullyPeriodic)
    ));
    assert_eq!(m, before);
}

#[test]
fn waters_are_whole_shared_classified_and_periodically_nonclashing() {
    let radius = Element::from_symbol("C")
        .unwrap()
        .van_der_waals_radius_angstrom(VanDerWaalsRadiusSource::Reference)
        .unwrap()
        * 0.1;
    let water_radius = 0.315_075_240_657_512_4 * 0.5612310241546865;
    for shape in [BoxShape::Cube, BoxShape::Dodecahedron, BoxShape::Octahedron] {
        let mut m = model("C");
        m.add_periodic_box(&PeriodicBoxOptions::Padding {
            padding: Quantity::new(1.0, NANOMETER),
            shape,
        })
        .unwrap();
        let old = m.clone();
        let report = m.add_solvent(&options()).unwrap();
        assert!(report.waters_added > 100);
        assert_eq!(m.atom_count(), 1 + 3 * report.waters_added);
        assert_eq!(m.bond_count(), 2 * report.waters_added);
        assert_eq!(m.topology().definition_count(), 2);
        assert_eq!(m.cell(), old.cell());
        assert_eq!(points(&m)[0], points(&old)[0]);
        assert!(!Arc::ptr_eq(&m.shared_topology(), &old.shared_topology()));
        let geometry = PeriodicGeometry::new(*m.cell().unwrap()).unwrap();
        let mut oxygens = vec![];
        for water in m.topology().molecules().skip(1) {
            assert_eq!(water.class(), MoleculeClass::Water);
            assert_eq!(water.molecule().atom_count(), 3);
            assert_eq!(water.molecule().bond_count(), 2);
            let atoms: Vec<_> = water.atoms().collect();
            assert_eq!(atoms[0].1.element.symbol(), "O");
            let p: Vec<_> = atoms
                .iter()
                .map(|(id, _)| m.position(*id).unwrap().into_value())
                .collect();
            for h in &p[1..] {
                assert!(((*h - p[0]).norm() - 0.09572).abs() < 0.0002);
            }
            let cos =
                (p[1] - p[0]).dot(p[2] - p[0]) / ((p[1] - p[0]).norm() * (p[2] - p[0]).norm());
            assert!((cos.acos().to_degrees() - 104.52).abs() < 0.2);
            assert!(
                geometry
                    .minimum_image(p[0] - Point3::origin())
                    .unwrap()
                    .norm()
                    >= radius + water_radius - 1e-12
            );
            for previous in &oxygens {
                assert!(
                    geometry.minimum_image(p[0] - *previous).unwrap().norm()
                        >= water_radius - 1e-12
                );
            }
            oxygens.push(p[0]);
        }
        for residue in m.residues() {
            assert_eq!(residue.class(), ResidueClass::Water);
        }
        assert_eq!(m.atom_sites().count(), 3 * report.waters_added);
    }
}

#[test]
fn charge_salt_rounding_and_seeded_ion_placement() {
    for (source, charge) in [("[Na+]", 1_i64), ("[Cl-]", -1), ("[Mg+2]", 2), ("C", 0)] {
        let initial = boxed(source);
        let mut water_only = initial.clone();
        let nw = water_only.add_solvent(&options()).unwrap().waters_added;
        let opts = SolventOptions {
            ionic_strength: Quantity::new(0.2, MOLAR),
            ..Default::default()
        };
        let mut a = initial.clone();
        let mut b = initial.clone();
        let report = a.add_solvent(&opts).unwrap();
        assert_eq!(b.add_solvent(&opts).unwrap(), report);
        assert_eq!(a.positions(), b.positions());
        let pairs =
            (((nw - charge.unsigned_abs() as usize) as f64 * 0.2 / 55.4) + 0.5).floor() as usize;
        assert_eq!(
            report.positive_ions_added,
            pairs + (-charge).max(0) as usize
        );
        assert_eq!(report.negative_ions_added, pairs + charge.max(0) as usize);
        assert_eq!(
            report.waters_added + report.positive_ions_added + report.negative_ions_added,
            nw
        );
        assert_eq!(
            a.atoms()
                .map(|(_, atom)| i64::from(atom.formal_charge))
                .sum::<i64>(),
            0
        );
        let geom = PeriodicGeometry::new(*a.cell().unwrap()).unwrap();
        let ions: Vec<_> = a
            .topology()
            .molecules()
            .filter(|m| m.class() == MoleculeClass::Ion)
            .map(|m| {
                a.position(m.atoms().next().unwrap().0)
                    .unwrap()
                    .into_value()
            })
            .collect();
        for i in 0..ions.len() {
            for j in i + 1..ions.len() {
                assert!(geom.minimum_image(ions[i] - ions[j]).unwrap().norm() > 0.5);
            }
        }
        let mut changed = initial;
        changed
            .add_solvent(&SolventOptions { seed: 42, ..opts })
            .unwrap();
        assert_ne!(changed.positions(), a.positions());
    }
}

#[test]
fn overrides_validate_topology_units_charge_and_ion_species() {
    let mut m = boxed("C");
    let radii =
        SolvationRadii::new(m.shared_topology(), Quantity::new(vec![2.0], ANGSTROM)).unwrap();
    let other = boxed("C");
    let wrong =
        SolvationRadii::new(other.shared_topology(), Quantity::new(vec![0.2], NANOMETER)).unwrap();
    let before = m.clone();
    assert!(matches!(
        m.add_solvent(&SolventOptions {
            solute_radii: Some(wrong),
            ..Default::default()
        }),
        Err(SolvationError::RadiusTopologyMismatch)
    ));
    assert_eq!(m, before);
    assert!(SolvationRadii::new(m.shared_topology(), Quantity::new(vec![], NANOMETER)).is_err());
    for r in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(
            SolvationRadii::new(m.shared_topology(), Quantity::new(vec![r], NANOMETER)).is_err()
        );
    }
    assert!(SolvationRadii::new(m.shared_topology(), Quantity::new(vec![1.0], MOLAR)).is_err());
    for (pos, neg) in [
        (PositiveIon::Lithium, NegativeIon::Fluoride),
        (PositiveIon::Potassium, NegativeIon::Bromide),
        (PositiveIon::Rubidium, NegativeIon::Iodide),
        (PositiveIon::Cesium, NegativeIon::Chloride),
    ] {
        let mut m = before.clone();
        let report = m
            .add_solvent(&SolventOptions {
                solute_radii: Some(radii.clone()),
                net_charge: Some(-2),
                positive_ion: pos,
                negative_ion: neg,
                ionic_strength: Quantity::new(0.1, MOLAR),
                ..Default::default()
            })
            .unwrap();
        assert_eq!(report.positive_ions_added, report.negative_ions_added + 2);
        assert_eq!(report.charge_basis, SolvationChargeBasis::Override);
        assert_eq!(m.atoms().next().unwrap().1.formal_charge, 0);
    }
    let composed = MOLE.try_div(METER.try_powi(3).unwrap()).unwrap();
    assert_eq!(
        Quantity::new(1000.0, composed)
            .into_unit(MOLAR)
            .unwrap()
            .into_value(),
        1.0
    );
}

#[test]
fn errors_leave_every_part_of_model_unchanged() {
    let mut no_cell = model("C");
    let before = no_cell.clone();
    assert!(matches!(
        no_cell.add_solvent(&options()),
        Err(SolvationError::MissingCell)
    ));
    assert_eq!(no_cell, before);
    let mut m = boxed("C");
    let before = m.clone();
    let bad = [
        SolventOptions {
            ionic_strength: Quantity::new(-1.0, MOLAR),
            ..options()
        },
        SolventOptions {
            ionic_strength: Quantity::new(f64::NAN, MOLAR),
            ..options()
        },
        SolventOptions {
            ionic_strength: Quantity::new(0.1, NANOMETER),
            ..options()
        },
        SolventOptions {
            ion_separation: Quantity::new(0.0, NANOMETER),
            ..options()
        },
        SolventOptions {
            max_candidates: 1,
            ..options()
        },
        SolventOptions {
            net_charge: Some(i64::MAX),
            ..Default::default()
        },
        SolventOptions {
            net_charge: Some(1),
            ion_separation: Quantity::new(100.0, NANOMETER),
            ..Default::default()
        },
        SolventOptions {
            ionic_strength: Quantity::new(1000.0, MOLAR),
            ..options()
        },
    ];
    for opts in bad {
        assert!(m.add_solvent(&opts).is_err());
        assert_eq!(m, before);
    }
    m.set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(3.0, 3.0, 3.0), NANOMETER),
            [true, true, false],
        )
        .unwrap(),
    ));
    let before = m.clone();
    assert!(matches!(
        m.add_solvent(&options()),
        Err(SolvationError::NotFullyPeriodic)
    ));
    assert_eq!(m, before);
    m.set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(0.1, 3.0, 3.0), NANOMETER),
            [true; 3],
        )
        .unwrap(),
    ));
    let before = m.clone();
    assert!(matches!(
        m.add_solvent(&options()),
        Err(SolvationError::NoSolventSpace)
    ));
    assert_eq!(m, before);
}

#[test]
fn original_entities_properties_and_existing_water_survive() {
    let mut builder = model("[H]O[H]").into_builder();
    let ids: Vec<_> = builder.topology_builder().atom_ids().collect();
    let chain = builder.hierarchy_mut().add_chain("SOL1", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "HOH", None, Some("123".into()), None)
        .unwrap();
    for (i, id) in ids.iter().enumerate() {
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                *id,
                AtomSiteMetadata {
                    label_atom_id: Some(format!("old{i}")),
                    ..Default::default()
                },
            )
            .unwrap();
    }
    builder
        .topology_builder_mut()
        .insert_property(key("topology_note"), PropertyValue::Int(1))
        .unwrap();
    let mut m = builder.build().unwrap();
    m.insert_property(key("model_note"), PropertyValue::Int(2))
        .unwrap();
    m.set_atom_property(ids[0], key("atom_note"), Some(PropertyValue::Int(3)))
        .unwrap();
    m.add_periodic_box(&PeriodicBoxOptions::Dimensions(Quantity::new(
        Vector3::new(3.0, 3.0, 3.0),
        NANOMETER,
    )))
    .unwrap();
    let before = m.clone();
    let report = m.add_solvent(&options()).unwrap();
    assert_eq!(&m.topology().atom_ids()[..ids.len()], ids);
    assert_eq!(&points(&m)[..ids.len()], points(&before));
    assert_eq!(
        m.atom_property(ids[0], &key("atom_note")).unwrap(),
        Some(PropertyValue::Int(3))
    );
    assert_eq!(
        m.atom_property(m.topology().atom_ids()[3], &key("atom_note"))
            .unwrap(),
        None
    );
    assert_eq!(
        report.cleared_topology_properties,
        vec![key("topology_note")]
    );
    assert_eq!(report.cleared_model_properties, vec![key("model_note")]);
    assert_eq!(
        m.hierarchy().chain(chain).unwrap(),
        before.hierarchy().chain(chain).unwrap()
    );
    assert_eq!(
        m.hierarchy().residue(residue).unwrap(),
        before.hierarchy().residue(residue).unwrap()
    );
    assert!(m.chains().any(|c| c.label_id() == "SOL2"));
    assert_eq!(before.topology().instance_count(), 1);
}

#[test]
fn externally_supplied_solute_and_increasing_boxes() {
    // Existing externally sourced CIP corpus molecule; no generated benchmark solute.
    let doc = kekule::sdf::parse_str(include_str!("fixtures/cip/VS132.sdf")).unwrap();
    let solute = doc.records()[0].to_model().unwrap();
    let mut previous = 0;
    for size in [3.0, 4.0, 6.0] {
        let mut m = solute.clone();
        m.add_periodic_box(&PeriodicBoxOptions::Dimensions(Quantity::new(
            Vector3::new(size, size, size),
            NANOMETER,
        )))
        .unwrap();
        let report = m.add_solvent(&options()).unwrap();
        assert!(report.waters_added > previous);
        let density = report.waters_added as f64 / size.powi(3);
        assert!((25.0..35.0).contains(&density), "density={density}");
        assert_eq!(&points(&m)[..solute.atom_count()], points(&solute));
        previous = report.waters_added;
    }
}

#[test]
fn rotated_skewed_and_left_handed_cells_fill_without_moving_solute() {
    for vectors in [
        [
            Vector3::new(0.0, 2.0, 0.0),
            Vector3::new(2.5, 0.0, 0.0),
            Vector3::new(0.2, 0.4, 2.0),
        ],
        [
            Vector3::new(2.0, 0.3, 0.1),
            Vector3::new(1.8, 1.5, 0.2),
            Vector3::new(0.3, 0.4, 2.0),
        ],
    ] {
        let mut m = model("C");
        let atom = m.topology().atom_ids()[0];
        let position = Point3::new(25.0, -11.0, 8.0);
        m.set_atom_positions([(atom, Quantity::new(position, NANOMETER))])
            .unwrap();
        let cell = PeriodicCell::new(Quantity::new(vectors, NANOMETER), [true; 3]).unwrap();
        m.add_periodic_box(&PeriodicBoxOptions::Cell(cell)).unwrap();
        let report = m.add_solvent(&options()).unwrap();
        assert!(report.waters_added > 100);
        assert_eq!(m.position(atom).unwrap().into_value(), position);
        let geometry = PeriodicGeometry::new(cell).unwrap();
        for water in m.topology().molecules().skip(1) {
            let oxygen = water.atoms().next().unwrap().0;
            let f = geometry
                .fractional(m.position(oxygen).unwrap().into_value() - position)
                .unwrap();
            assert!(f.iter().all(|v| *v >= -0.5 - 1e-12 && *v < 0.5 + 1e-12));
        }
    }
}

#[test]
fn core_periodic_kernel_checks_limits_and_shortest_skew_translation() {
    let cell = PeriodicCell::new(
        Quantity::new(
            [
                Vector3::new(2.0, 0.0, 0.0),
                Vector3::new(1.9, 0.3, 0.0),
                Vector3::new(0.0, 0.0, 2.0),
            ],
            NANOMETER,
        ),
        [true; 3],
    )
    .unwrap();
    let geometry = PeriodicGeometry::new(cell).unwrap();
    assert!((geometry.shortest_translation().unwrap().norm() - 0.1_f64.hypot(0.3)).abs() < 1e-12);
    assert!(matches!(
        geometry.minimum_image(Vector3::new(f64::NAN, 0.0, 0.0)),
        Err(PeriodicGeometryError::NumericalFailure)
    ));
    let bad = PeriodicCell::new(
        Quantity::new(
            [
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(1.0, 1e-10, 0.0),
                Vector3::new(0.0, 0.0, 1.0),
            ],
            NANOMETER,
        ),
        [true; 3],
    )
    .unwrap();
    assert!(matches!(
        PeriodicGeometry::new(bad).unwrap().shortest_translation(),
        Err(PeriodicGeometryError::ImageSearchLimit)
    ));
}

#[test]
fn absent_elemental_radius_requires_an_explicit_override() {
    let mut m = boxed("[Lv]");
    let before = m.clone();
    assert!(matches!(
        m.add_solvent(&options()),
        Err(SolvationError::MissingElementRadius(_))
    ));
    assert_eq!(m, before);
    let radii =
        SolvationRadii::new(m.shared_topology(), Quantity::new(vec![0.2], NANOMETER)).unwrap();
    m.add_solvent(&SolventOptions {
        solute_radii: Some(radii.clone()),
        ..options()
    })
    .unwrap();
    let after = m.clone();
    assert!(matches!(
        m.add_solvent(&SolventOptions {
            solute_radii: Some(radii),
            ..options()
        }),
        Err(SolvationError::RadiusTopologyMismatch)
    ));
    assert_eq!(m, after);
}

#[test]
fn all_candidate_waters_can_be_replaced_without_unused_definitions() {
    let initial = boxed("C");
    let mut water_only = initial.clone();
    let count = water_only.add_solvent(&options()).unwrap().waters_added;
    let mut m = initial;
    let report = m
        .add_solvent(&SolventOptions {
            net_charge: Some(count as i64),
            ion_separation: Quantity::new(0.01, NANOMETER),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(report.waters_added, 0);
    assert_eq!(report.negative_ions_added, count);
    assert_eq!(m.topology().definition_count(), 2);
    assert_eq!(m.atom_count(), 1 + count);
    assert_eq!(m.bond_count(), 0);
}

#[test]
fn translation_that_loses_water_geometry_is_rejected_atomically() {
    let mut m = boxed("C");
    let id = m.topology().atom_ids()[0];
    m.set_atom_positions([(id, Quantity::new(Point3::new(1e20, 0.0, 0.0), NANOMETER))])
        .unwrap();
    let before = m.clone();
    assert!(matches!(
        m.add_solvent(&options()),
        Err(SolvationError::Geometry(
            PeriodicGeometryError::NumericalFailure
        ))
    ));
    assert_eq!(m, before);
}

#[test]
fn salt_rounds_half_up_and_molarity_accepts_composed_units() {
    let initial = boxed("C");
    let mut water_only = initial.clone();
    let count = water_only.add_solvent(&options()).unwrap().waters_added;
    for (fraction, expected) in [(0.49, 0), (0.5, 1), (1.49, 1), (1.5, 2)] {
        let concentration = fraction * 55.4 / count as f64;
        let mut m = initial.clone();
        let report = m
            .add_solvent(&SolventOptions {
                ionic_strength: Quantity::new(concentration, MOLAR),
                ..options()
            })
            .unwrap();
        assert_eq!(report.positive_ions_added, expected);
        assert_eq!(report.negative_ions_added, expected);
    }
    let composed = MOLE.try_div(METER.try_powi(3).unwrap()).unwrap();
    let mut a = initial.clone();
    let mut b = initial;
    let a_report = a
        .add_solvent(&SolventOptions {
            ionic_strength: Quantity::new(150.0, composed),
            ..options()
        })
        .unwrap();
    let b_report = b
        .add_solvent(&SolventOptions {
            ionic_strength: Quantity::new(0.15, MOLAR),
            ..options()
        })
        .unwrap();
    assert_eq!(a_report, b_report);
    assert_eq!(a.positions(), b.positions());
}
