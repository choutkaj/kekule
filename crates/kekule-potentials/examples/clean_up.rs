//! Clean-up of an edited ligand: minimize a group with the rest of its
//! molecule fixed.
//!
//! ```text
//! cargo run --release -p kekule-potentials --example clean_up -- ligand.sdf [environment.cif] [--whole]
//! ```
//!
//! The first SDF record must contain explicit hydrogens with 3D coordinates.
//! The ligand is appended to the optional mmCIF environment to form a larger
//! model. The example then mimics an edit: it displaces the ligand's last three
//! heavy atoms (all of them in smaller ligands) and their hydrogens, the "new
//! group", by up to 0.3 Å. Clean-up takes the molecule of the target atoms out
//! of the model, parameterizes it with OpenFF Rosemary and the bundled Ash
//! charges, minimizes the targets with the rest of the molecule fixed, and
//! writes the coordinates back. `--whole` cleans up the whole ligand instead of
//! the group.

use std::{collections::BTreeSet, env, error::Error, fs, ops::ControlFlow, time::Instant};

use kekule::geometry::{Point3, Vector3};
use kekule::structure::Model;
use kekule::topology::{AtomSelection, TopologyAtomIndex};
use kekule::units::{Quantity, CANONICAL_LENGTH_UNIT, NANOMETER};
use kekule::{mmcif, sdf};
use kekule_openff::{ForceField, NaglModel};
use kekule_potentials::{
    minimize_scoped_with_observer, openff::OpenFfPotential, MinimizationScope, MinimizeOptions,
    Potential,
};

fn main() -> Result<(), Box<dyn Error>> {
    let mut paths = Vec::new();
    let mut whole = false;
    for arg in env::args().skip(1) {
        match arg.as_str() {
            "--whole" => whole = true,
            _ => paths.push(arg),
        }
    }
    let usage = "usage: clean_up <ligand.sdf> [environment.cif] [--whole]";
    let ligand_path = paths.first().ok_or(usage)?;

    let document = sdf::parse_str(&fs::read_to_string(ligand_path)?, Default::default())?;
    let record = document.records().first().ok_or("no SDF records")?;
    let ligand = record.interpret(Default::default())?.into_model();
    let mut system = match paths.get(1) {
        Some(path) => {
            let environment = mmcif::parse_str(&fs::read_to_string(path)?, Default::default())?
                .interpret(Default::default())?
                .into_model();
            let mut editor = environment.to_editor();
            editor.append_model(&ligand)?;
            editor.finish()?
        }
        None => ligand.clone(),
    };
    // The ligand's instances were appended last.
    let first_ligand = system.topology().instance_count() - ligand.topology().instance_count();
    let ligand_atoms = AtomSelection::from_predicate(system.shared_topology(), |atom| {
        atom.id().instance().index() >= first_ligand
    });

    // Mimic a builder edit of the last three heavy atoms and their hydrogens.
    let heavy: Vec<_> = ligand_atoms
        .atom_ids()
        .filter(|&id| {
            system
                .topology()
                .atom(id)
                .unwrap()
                .element()
                .atomic_number()
                != 1
        })
        .collect();
    // Ligands with fewer than three heavy atoms edit all of them.
    let group: BTreeSet<_> = heavy[heavy.len().saturating_sub(3)..]
        .iter()
        .flat_map(|&id| {
            let atom = system.topology().atom(id).unwrap();
            std::iter::once(id).chain(
                atom.neighbors()
                    .filter(|n| n.element().atomic_number() == 1)
                    .map(|n| n.id()),
            )
        })
        .collect();
    let targets = if whole {
        ligand_atoms.clone()
    } else if group.is_empty() {
        return Err("the ligand has no heavy atoms to edit; pass --whole".into());
    } else {
        AtomSelection::from_atoms(system.shared_topology(), group.iter().copied())?
    };
    for (n, id) in group.iter().enumerate() {
        let p = system.position(*id).unwrap().value_in(NANOMETER)?;
        let k = n as f64 + 1.0;
        let shift = Vector3::new((1.3 * k).sin(), (2.1 * k).cos(), (0.7 * k).sin()) * 0.017;
        system.set_position(*id, Quantity::new(p + shift, NANOMETER))?;
    }
    let edited = system.clone();
    println!(
        "system: {} atoms in {} molecules; ligand {} atoms; {} target atoms",
        system.atom_count(),
        system.topology().instance_count(),
        ligand_atoms.len(),
        targets.len()
    );

    // Clean-up: the molecules of the targets, the targets moving.
    let clock = Instant::now();
    let molecules = targets.expand_to_instances();
    let molecule: Model = system.subset(&molecules)?;
    // Subsets keep dense order, so subset atom `k` is the `k`-th selected atom.
    let moving = molecules
        .indices()
        .iter()
        .enumerate()
        .filter(|(_, index)| targets.contains_index(**index))
        .map(|(k, _)| TopologyAtomIndex::new(k as u32));
    let scope = MinimizationScope::moving(AtomSelection::from_indices(
        molecule.shared_topology(),
        moving,
    )?);
    let subset_time = clock.elapsed();

    // Loading is a one-off cost; keep the force field and charge model.
    let clock = Instant::now();
    let force_field = ForceField::rosemary()?;
    let nagl = NaglModel::ash()?;
    let load_time = clock.elapsed();

    let clock = Instant::now();
    let parameters = force_field.parameterize(molecule.shared_topology(), &nagl)?;
    let potential = OpenFfPotential::new(&parameters)?;
    let parameter_time = clock.elapsed();
    let charge: f64 = parameters.charges().value().iter().sum();
    println!("net partial charge: {charge:+.6} e");

    check_gradient(&potential, &molecule)?;

    let clock = Instant::now();
    let mut steps = 0;
    let result = minimize_scoped_with_observer(
        &potential,
        molecule.as_model_view(),
        &scope,
        // Whole drug-like ligands from ideal geometries can take over 1000 steps.
        &MinimizeOptions::default().with_max_iterations(5_000),
        |_| {
            steps += 1;
            // A GUI returns `Break` here when the user cancels.
            ControlFlow::Continue(())
        },
    )?;
    let minimize_time = clock.elapsed();

    // Write the moving atoms back in one edit.
    let minimized = result.positions().values().into_value();
    for (k, index) in molecules.indices().iter().enumerate() {
        if targets.contains_index(*index) {
            let id = system.topology().atom_ids()[index.index()];
            system.set_position(id, Quantity::new(minimized[k], CANONICAL_LENGTH_UNIT))?;
        }
    }
    let unchanged = system
        .atoms()
        .zip(edited.atoms())
        .filter(|(after, _)| !targets.contains(after.id()))
        .all(|(after, before)| bits(after.position()) == bits(before.position()));

    println!(
        "potential energy: {:.3} -> {:.3} kJ/mol ({:?}, {} steps, {} evaluations, max moving gradient {:.3} kJ/(mol nm))",
        result.initial_energy().total().into_value(),
        result.final_evaluation().energy().total().into_value(),
        result.status(),
        result.iterations(),
        result.evaluations(),
        result.max_gradient().into_value(),
    );
    assert_eq!(steps, result.iterations());
    println!("atoms outside the targets unchanged to the bit: {unchanged}");
    println!(
        "time: load {:.1} ms, subset {:.1} ms, parameterize {:.1} ms, minimize {:.1} ms",
        load_time.as_secs_f64() * 1e3,
        subset_time.as_secs_f64() * 1e3,
        parameter_time.as_secs_f64() * 1e3,
        minimize_time.as_secs_f64() * 1e3
    );
    Ok(())
}

fn bits(position: Quantity<Point3>) -> [u64; 3] {
    let p = position.into_value();
    [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
}

/// Compares the analytic gradient with central differences.
fn check_gradient(potential: &OpenFfPotential, model: &Model) -> Result<(), Box<dyn Error>> {
    let analytic = potential.evaluate(model.as_model_view())?;
    let h = 1e-6;
    let mut worst: f64 = 0.0;
    let mut largest: f64 = 0.0;
    let mut trial = model.clone();
    for (index, id) in model.topology().atom_ids().iter().enumerate() {
        let p = model.position(*id).unwrap().value_in(NANOMETER)?;
        let g = analytic.gradient().into_value()[index];
        for (axis, step) in [
            Vector3::new(h, 0.0, 0.0),
            Vector3::new(0.0, h, 0.0),
            Vector3::new(0.0, 0.0, h),
        ]
        .into_iter()
        .enumerate()
        {
            trial.set_position(*id, Quantity::new(p + step, NANOMETER))?;
            let plus = potential
                .energy(trial.as_model_view())?
                .total()
                .into_value();
            trial.set_position(*id, Quantity::new(p - step, NANOMETER))?;
            let minus = potential
                .energy(trial.as_model_view())?
                .total()
                .into_value();
            let numeric = (plus - minus) / (2.0 * h);
            let exact = [g.x, g.y, g.z][axis];
            worst = worst.max((numeric - exact).abs());
            largest = largest.max(exact.abs());
        }
        trial.set_position(*id, Quantity::new(p, NANOMETER))?;
    }
    println!(
        "gradient check: max |analytic - central difference| {worst:.2e} kJ/(mol nm), max |gradient| {largest:.1}"
    );
    Ok(())
}
