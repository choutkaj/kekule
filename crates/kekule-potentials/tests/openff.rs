#![cfg(feature = "openff")]

use std::sync::Arc;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::structure::{Ensemble, Model, ModelView, Positions};
use kekule::structure::{Trajectory, TrajectoryFrame};
use kekule::topology::{Topology, TopologyBuilder};
use kekule::units::{Quantity, ELEMENTARY_CHARGE, NANOMETER};
use kekule_openff::{ForceField, ParameterizedTopology};
use kekule_potentials::openff::{OpenFfPotential, OpenFfPotentialError};
use kekule_potentials::{
    minimize, ComponentKind, EvaluationError, MinimizationStatus, MinimizeOptions, Potential,
    SingularGeometry,
};
use kekule_traj::FrameBuffer;

// Deliberately simple, distinct numbers make every handler visible in the
// energy; the reference below recomputes it independently of the kernels.
const TOY: &str = r#"<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
<Bonds version="0.4"><Bond smirks="[*:1]~[*:2]" id="b" length="0.11*nanometer" k="200000*kilojoule_per_mole/nanometer**2"/></Bonds>
<Angles version="0.3"><Angle smirks="[*:1]~[*:2]~[*:3]" id="a" angle="110*degree" k="150*kilojoule_per_mole/radian**2"/></Angles>
<ProperTorsions version="0.4"><Proper smirks="[*:1]~[*:2]~[*:3]~[*:4]" id="t" periodicity1="3" phase1="0*degree" k1="1.5*kilojoule_per_mole" periodicity2="1" phase2="45*degree" k2="0.4*kilojoule_per_mole" idivf2="1"/></ProperTorsions>
<ImproperTorsions version="0.3"><Improper smirks="[*:1]~[#6X3:2](~[*:3])~[*:4]" id="i" periodicity1="2" phase1="180*degree" k1="4.5*kilojoule_per_mole"/></ImproperTorsions>
<vdW version="0.4"><Atom smirks="[*:1]" id="v" sigma="0.3*nanometer" epsilon="0.2*kilojoule_per_mole"/><Atom smirks="[#1:1]" id="vh" sigma="0.25*nanometer" epsilon="0.05*kilojoule_per_mole"/></vdW>
<Electrostatics version="0.4" scale14="0.8333333333"/>
<LibraryCharges version="0.3">
<LibraryCharge smirks="[#1:1]" id="qh" charge1="0.1*elementary_charge"/>
<LibraryCharge smirks="[#6X4:1]" id="qc4" charge1="-0.3*elementary_charge"/>
<LibraryCharge smirks="[#6X3:1]" id="qc3" charge1="0.3*elementary_charge"/>
<LibraryCharge smirks="[#8:1]" id="qo" charge1="-0.4*elementary_charge"/>
</LibraryCharges>
</SMIRNOFF>"#;

/// Acetaldehyde with explicit hydrogens: H0 C1 H2 H3 C4 H5 O6.
const ACETALDEHYDE: &str = "[H]C([H])([H])C([H])=O";
const GEOMETRY: [[f64; 3]; 7] = [
    [-0.036, 0.101, 0.002],
    [0.0, 0.0, 0.0],
    [-0.037, -0.052, 0.088],
    [-0.034, -0.049, -0.090],
    [0.151, 0.003, -0.004],
    [0.207, -0.092, 0.010],
    [0.214, 0.108, 0.001],
];

fn topology(instances: usize) -> Arc<Topology> {
    let mut molecule = kekule::smiles::to_molecules(ACETALDEHYDE)
        .unwrap()
        .remove(0);
    molecule.perceive().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(molecule.clone()).unwrap();
    for _ in 0..instances {
        builder.add_instance(definition).unwrap();
    }
    Arc::new(builder.build().unwrap())
}

/// Instances are offset along z so intermolecular pairs are evaluated at full strength.
fn model(topology: &Arc<Topology>) -> Model {
    let points = (0..topology.atom_count())
        .map(|i| {
            let [x, y, z] = GEOMETRY[i % 7];
            Point3::new(x, y + 0.01 * (i / 7) as f64, z + 0.42 * (i / 7) as f64)
        })
        .collect::<Vec<_>>();
    Model::new(
        Arc::clone(topology),
        Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap()
}

fn parameterize(topology: &Arc<Topology>) -> ParameterizedTopology {
    ForceField::from_offxml(TOY)
        .unwrap()
        .parameterize(
            Arc::clone(topology),
            kekule_openff::ChargeMethod::LibraryOnly,
        )
        .unwrap()
}

fn sub(a: Point3, b: Point3) -> [f64; 3] {
    [a.x - b.x, a.y - b.y, a.z - b.z]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Independent SMIRNOFF energy from the public parameterization, using
/// textbook formulas: acos angles, the IUPAC/OpenMM signed dihedral
/// `atan2(b2_hat . (n1 x n2), n1 . n2)`, and an explicit pair loop that looks
/// exceptions up by atom identity.
fn reference_energy(p: &ParameterizedTopology, model: &Model) -> [f64; 6] {
    let x = |atom| model.position(atom).unwrap().into_value();
    let mut e = [0.0; 6];
    for b in p.bonds() {
        let r = norm(sub(x(b.atoms[0]), x(b.atoms[1])));
        e[0] += 0.5 * b.parameter.k.value() * (r - b.parameter.length.value()).powi(2);
    }
    for a in p.angles() {
        let u = sub(x(a.atoms[0]), x(a.atoms[1]));
        let v = sub(x(a.atoms[2]), x(a.atoms[1]));
        let theta = (dot(u, v) / (norm(u) * norm(v))).acos();
        e[1] += 0.5 * a.parameter.k.value() * (theta - a.parameter.angle.value()).powi(2);
    }
    for (slot, torsions) in [(2, p.proper_torsions()), (3, p.improper_torsions())] {
        for t in torsions {
            let [a, b, c, d] = t.atoms.map(x);
            let (b1, b2, b3) = (sub(b, a), sub(c, b), sub(d, c));
            let (n1, n2) = (cross(b1, b2), cross(b2, b3));
            let phi = dot(cross(n1, n2), b2.map(|v| v / norm(b2))).atan2(dot(n1, n2));
            for term in &t.parameter.terms {
                e[slot] += term.k.value() / term.idivf
                    * (1.0 + (f64::from(term.periodicity) * phi - term.phase.value()).cos());
            }
        }
    }
    let ids = model.topology().atom_ids();
    let charges = p.charges().value();
    for i in 0..ids.len() {
        for j in i + 1..ids.len() {
            let exception = p
                .pair_exceptions()
                .iter()
                .find(|e| e.atoms == [ids[i], ids[j]] || e.atoms == [ids[j], ids[i]]);
            let (vs, qs) = exception.map_or((1.0, 1.0), |e| (e.vdw_scale, e.electrostatics_scale));
            let r = norm(sub(x(ids[i]), x(ids[j])));
            let (a, b) = (&p.vdw()[i], &p.vdw()[j]);
            let sigma = (a.sigma.value() + b.sigma.value()) / 2.0;
            let epsilon = (a.epsilon.value() * b.epsilon.value()).sqrt();
            e[4] += vs * 4.0 * epsilon * ((sigma / r).powi(12) - (sigma / r).powi(6));
            e[5] += qs * 138.935_457_644_381_98 * charges[i] * charges[j] / r;
        }
    }
    e
}

#[test]
fn energy_components_match_an_independent_smirnoff_reference() {
    let topology = topology(2);
    let p = parameterize(&topology);
    assert!(!p.improper_torsions().is_empty() && !p.pair_exceptions().is_empty());
    let model = model(&topology);
    let potential = OpenFfPotential::new(&p).unwrap();
    let energy = potential.energy(model.as_model_view()).unwrap();
    let expected = reference_energy(&p, &model);
    let kinds = [
        ComponentKind::Bonds,
        ComponentKind::Angles,
        ComponentKind::ProperTorsions,
        ComponentKind::ImproperTorsions,
        ComponentKind::VanDerWaals,
        ComponentKind::Electrostatics,
    ];
    assert_eq!(
        energy
            .components()
            .iter()
            .map(|c| c.kind)
            .collect::<Vec<_>>(),
        kinds
    );
    for (kind, expected) in kinds.into_iter().zip(expected) {
        let actual = energy.component(kind).unwrap().into_value();
        assert!(expected != 0.0, "{kind} must contribute");
        assert!(
            (actual - expected).abs() <= 1e-10 * (1.0 + expected.abs()),
            "{kind}: {actual} != {expected}"
        );
    }
    let sum: f64 = expected.iter().sum();
    assert!((energy.total().into_value() - sum).abs() <= 1e-10 * (1.0 + sum.abs()));

    let evaluation = potential.evaluate(model.as_model_view()).unwrap();
    assert_eq!(evaluation.energy(), &energy);
    assert!(Arc::ptr_eq(evaluation.topology(), &topology));
}

#[test]
fn gradient_differentiates_energy_and_components_add_up() {
    let topology = topology(2);
    let potential = OpenFfPotential::new(&parameterize(&topology)).unwrap();
    let model = model(&topology);
    let evaluation = potential.evaluate(model.as_model_view()).unwrap();
    let gradient = evaluation.gradient().into_value();
    let points = model.positions().values().into_value().to_vec();
    let energy_at = |points: &[Point3]| {
        let mut moved = model.clone();
        moved
            .conformation_mut()
            .set_positions(Quantity::new(points, NANOMETER))
            .unwrap();
        potential
            .energy(moved.as_model_view())
            .unwrap()
            .total()
            .into_value()
    };
    for atom in 0..points.len() {
        for axis in 0..3 {
            let estimate = |h: f64| {
                let shift = |sign: f64| {
                    let mut moved = points.clone();
                    let p = &mut moved[atom];
                    *[&mut p.x, &mut p.y, &mut p.z][axis] += sign * h;
                    energy_at(&moved)
                };
                (shift(1.0) - shift(-1.0)) / (2.0 * h)
            };
            let numerical = (4.0 * estimate(5e-6) - estimate(1e-5)) / 3.0;
            let g = gradient[atom];
            let analytic = [g.x, g.y, g.z][axis];
            assert!(
                (analytic - numerical).abs() <= 1e-5 * (1.0 + numerical.abs()),
                "atom {atom} axis {axis}: {analytic} != {numerical}"
            );
        }
    }

    let components = potential
        .evaluate_components(model.as_model_view())
        .unwrap();
    for (atom, total) in gradient.iter().enumerate() {
        let mut sum = Vector3::zero();
        for c in &components {
            let g = c.gradient.value()[atom];
            sum = Vector3::new(sum.x + g.x, sum.y + g.y, sum.z + g.z);
        }
        assert!(
            norm(sub(
                Point3::new(sum.x, sum.y, sum.z),
                Point3::new(total.x, total.y, total.z)
            )) < 1e-9
        );
    }
    for (component, reported) in components.iter().zip(evaluation.energy().components()) {
        assert_eq!(component.kind, reported.kind);
        assert_eq!(component.energy, reported.energy);
    }
}

#[test]
fn energy_is_rigid_motion_invariant_and_gradient_has_no_net_force() {
    let topology = topology(2);
    let potential = OpenFfPotential::new(&parameterize(&topology)).unwrap();
    let model = model(&topology);
    let reference = potential.evaluate(model.as_model_view()).unwrap();
    let net = reference
        .gradient()
        .into_value()
        .iter()
        .fold([0.0; 3], |n, g| [n[0] + g.x, n[1] + g.y, n[2] + g.z]);
    assert!(norm(net) < 1e-9, "net gradient {net:?}");

    let (c, s) = (0.3_f64.cos(), 0.3_f64.sin());
    let moved = model
        .positions()
        .values()
        .into_value()
        .iter()
        .map(|p| Point3::new(c * p.x - s * p.y + 1.7, s * p.x + c * p.y - 0.4, p.z + 2.2))
        .collect::<Vec<_>>();
    let mut transformed = model.clone();
    transformed
        .conformation_mut()
        .set_positions(Quantity::new(moved, NANOMETER))
        .unwrap();
    let total = |m: &Model| {
        potential
            .energy(m.as_model_view())
            .unwrap()
            .total()
            .into_value()
    };
    assert!((total(&transformed) - total(&model)).abs() < 1e-10 * total(&model).abs().max(1.0));
}

#[test]
fn potential_binds_its_layout_and_rejects_periodic_cells() {
    let topology = topology(1);
    let potential = OpenFfPotential::new(&parameterize(&topology)).unwrap();
    assert!(Arc::ptr_eq(potential.topology(), &topology));
    let model = model(&topology);
    // An independently published equal topology is a different layout.
    let other = self::model(&self::topology(1));
    assert_eq!(
        potential.energy(other.as_model_view()),
        Err(EvaluationError::IncompatibleTopology)
    );
    // Reperception publishes a new snapshot of the same layout, so the
    // prepared potential evaluates both identically.
    let mut perceived = model.clone();
    perceived.perceive().unwrap();
    let evaluation = potential.evaluate(perceived.as_model_view()).unwrap();
    assert!(Arc::ptr_eq(
        evaluation.topology(),
        &perceived.shared_topology()
    ));
    assert_eq!(
        evaluation.energy().total(),
        potential.energy(model.as_model_view()).unwrap().total()
    );

    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(3.0, 3.0, 3.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    let mut periodic = model.clone();
    periodic.conformation_mut().set_cell(Some(cell));
    assert_eq!(
        potential.energy(periodic.as_model_view()),
        Err(EvaluationError::UnsupportedPeriodicCell)
    );
    let ensemble = Ensemble::from_models([periodic]).unwrap();
    assert_eq!(
        potential.evaluate(ensemble.get(0).unwrap().as_model_view()),
        Err(EvaluationError::UnsupportedPeriodicCell)
    );
}

#[test]
fn ensemble_members_and_trajectory_frames_evaluate_like_models() {
    let topology = topology(1);
    let potential = OpenFfPotential::new(&parameterize(&topology)).unwrap();
    let model = model(&topology);
    let expected = potential.evaluate(model.as_model_view()).unwrap();

    let ensemble = Ensemble::from_models([model.clone()]).unwrap();
    assert_eq!(
        potential
            .evaluate(ensemble.get(0).unwrap().as_model_view())
            .unwrap(),
        expected
    );
    let frame = TrajectoryFrame::new(model.positions().clone());
    let trajectory = Trajectory::from_items(Arc::clone(&topology), [frame]).unwrap();
    assert_eq!(
        potential
            .evaluate(trajectory.get(0).unwrap().as_model_view())
            .unwrap(),
        expected
    );
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    buffer
        .frame_mut()
        .conformation_mut()
        .set_positions(model.positions().values())
        .unwrap();
    assert_eq!(
        potential.evaluate(buffer.as_model_view()).unwrap(),
        expected
    );
}

#[test]
fn singular_coordinates_name_the_qualified_atoms() {
    let topology = topology(2);
    let potential = OpenFfPotential::new(&parameterize(&topology)).unwrap();
    let mut model = model(&topology);
    let ids = model.topology().atom_ids().to_vec();
    // Place one instance's oxygen onto the other instance's methyl carbon.
    let carbon = model.position(ids[1]).unwrap();
    model.set_position(ids[13], carbon).unwrap();
    let expected = EvaluationError::InvalidGeometry {
        interaction: "nonbonded pair",
        atoms: vec![ids[1], ids[13]],
        kind: SingularGeometry::CoincidentAtoms,
    };
    assert_eq!(
        potential.energy(model.as_model_view()),
        Err(expected.clone())
    );
    assert_eq!(potential.evaluate(model.as_model_view()), Err(expected));
}

#[test]
fn replacement_charges_change_only_electrostatics() {
    let topology = topology(1);
    let p = parameterize(&topology);
    let model = model(&topology);
    let assigned = OpenFfPotential::new(&p).unwrap();
    let neutral = OpenFfPotential::new(&p)
        .unwrap()
        .with_charges(Quantity::new(vec![0.0; 7], ELEMENTARY_CHARGE))
        .unwrap();
    assert_eq!(neutral.charges().value(), &[0.0; 7]);
    let (a, n) = (
        assigned.energy(model.as_model_view()).unwrap(),
        neutral.energy(model.as_model_view()).unwrap(),
    );
    for c in a.components() {
        let other = n.component(c.kind).unwrap();
        if c.kind == ComponentKind::Electrostatics {
            assert_ne!(c.energy, other);
            assert_eq!(other.into_value(), 0.0);
        } else {
            assert_eq!(c.energy, other);
        }
    }
    assert!(matches!(
        OpenFfPotential::new(&p)
            .unwrap()
            .with_charges(Quantity::new(vec![0.0; 6], ELEMENTARY_CHARGE)),
        Err(OpenFfPotentialError::AtomCountMismatch {
            expected: 7,
            actual: 6,
            ..
        })
    ));
    assert!(OpenFfPotential::new(&p)
        .unwrap()
        .with_charges(Quantity::new(vec![0.0; 7], NANOMETER))
        .is_err());
}

#[test]
fn nonperiodic_vdw_cutoff_is_rejected_at_preparation() {
    let xml = TOY.replace(
        r#"<vdW version="0.4">"#,
        r#"<vdW version="0.4" nonperiodic_method="cutoff">"#,
    );
    let topology = topology(1);
    let p = ForceField::from_offxml(&xml)
        .unwrap()
        .parameterize(
            Arc::clone(&topology),
            kekule_openff::ChargeMethod::LibraryOnly,
        )
        .unwrap();
    assert_eq!(
        OpenFfPotential::new(&p).unwrap_err(),
        OpenFfPotentialError::UnsupportedNonbondedMethod {
            handler: "vdW",
            method: "cutoff",
        }
    );
}

#[test]
fn minimization_reaches_a_stationary_point_without_changing_its_input() {
    let topology = topology(2);
    let potential = OpenFfPotential::new(&parameterize(&topology)).unwrap();
    let model = model(&topology);
    let original = model.clone();
    let options = MinimizeOptions {
        gradient_tolerance: Quantity::new(1e-3, kekule::units::KILOJOULE_PER_MOLE_PER_NANOMETER),
        max_iterations: 5_000,
        ..MinimizeOptions::default()
    };
    let result = minimize(&potential, model.as_model_view(), &options).unwrap();
    assert_eq!(result.status(), MinimizationStatus::Converged);
    assert_eq!(model, original);
    let initial = result.initial_energy().total().into_value();
    let last = result.final_evaluation();
    assert!(last.energy().total().into_value() < initial);
    assert!(last.max_gradient_norm().into_value() <= 1e-3);

    let minimized = result.to_model(model.as_model_view()).unwrap();
    assert!(Arc::ptr_eq(&minimized.shared_topology(), &topology));
    assert_eq!(minimized.positions(), result.positions());
    assert_eq!(
        &potential.evaluate(minimized.as_model_view()).unwrap(),
        result.final_evaluation()
    );
    // Bonds relax to their 0.11 nm equilibrium within a small strain.
    for b in parameterize(&topology).bonds() {
        let r = norm(sub(
            minimized.position(b.atoms[0]).unwrap().into_value(),
            minimized.position(b.atoms[1]).unwrap().into_value(),
        ));
        assert!((r - 0.11).abs() < 0.01, "bond length {r}");
    }
    let unrelated = self::model(&self::topology(2));
    assert!(result.to_model(unrelated.as_model_view()).is_err());
}

#[test]
fn potentials_are_shareable_trait_objects() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<OpenFfPotential>();
    let topology = topology(1);
    let potential: Box<dyn Potential> =
        Box::new(OpenFfPotential::new(&parameterize(&topology)).unwrap());
    let model = model(&topology);
    let view: ModelView<'_> = model.as_model_view();
    let result = minimize(potential.as_ref(), view, &MinimizeOptions::default()).unwrap();
    assert!(result.iterations() > 0);
}

#[test]
fn atom_parameters_follow_an_interleaved_dense_order() {
    use kekule::core::AtomId;
    use kekule::topology::{InstanceAtomId, MoleculeInstanceId};

    let ordered = topology(2);
    // Alternate the two instances atom by atom, as a source file's rows may.
    let interleaved_order = (0..7).flat_map(|local| {
        [0, 1].map(|instance| {
            InstanceAtomId::new(MoleculeInstanceId::new(instance), AtomId::new(local))
        })
    });
    let mut builder = Arc::try_unwrap(topology(2)).unwrap().into_builder();
    builder.set_atom_order(interleaved_order).unwrap();
    let interleaved = Arc::new(builder.build().unwrap());
    let reference = model(&ordered);
    let moved = Model::new(
        Arc::clone(&interleaved),
        Positions::new(Quantity::new(
            interleaved
                .atom_ids()
                .iter()
                .map(|&atom| reference.position(atom).unwrap().into_value())
                .collect::<Vec<_>>(),
            NANOMETER,
        ))
        .unwrap(),
    )
    .unwrap();

    let (source, target) = (parameterize(&ordered), parameterize(&interleaved));
    for (dense, &atom) in interleaved.atom_ids().iter().enumerate() {
        let original = ordered.atom_index(atom).unwrap().index();
        assert_eq!(target.vdw()[dense].sigma, source.vdw()[original].sigma);
        assert_eq!(target.vdw()[dense].epsilon, source.vdw()[original].epsilon);
        assert_eq!(
            target.charges().value()[dense],
            source.charges().value()[original]
        );
    }
    let expected = OpenFfPotential::new(&source)
        .unwrap()
        .energy(reference.as_model_view())
        .unwrap()
        .total()
        .into_value();
    let actual = OpenFfPotential::new(&target)
        .unwrap()
        .energy(moved.as_model_view())
        .unwrap()
        .total()
        .into_value();
    assert!((actual - expected).abs() <= 1e-10 * (1.0 + expected.abs()));
}
