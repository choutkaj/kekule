//! Informational cost of non-covalent interaction searches in a structure of
//! about 100,000 atoms: 66 copies of KRas G12C with its inhibitor (PDB
//! 9UOH) packed 45 Å apart, searched from one copy's inhibitor, from it and
//! its pocket, and from a whole copy, within the structure and against the
//! inhibitor as a separate structure.

use std::error::Error;
use std::hint::black_box;
use std::sync::Arc;
use std::time::{Duration, Instant};

use kekule::geometry::{Matrix3, RigidTransform, Vector3};
use kekule::interactions::{
    find_interactions_between, find_interactions_within, InteractionCriteria, InteractionFeatures,
    InteractionSide, Interactions,
};
use kekule::topology::AtomSelection;
use kekule::units::{Quantity, ANGSTROM};

const COPIES: u32 = 66;
const SPACING: f64 = 45.0;
const PASSES: usize = 21;

fn main() -> Result<(), Box<dyn Error>> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/mmcif/9UOH.cif");
    let kras = kekule::mmcif::parse_str(&std::fs::read_to_string(path)?, Default::default())?
        .interpret(Default::default())?
        .into_model();
    let mut editor = kras.to_editor();
    for copy in 1..COPIES {
        let offset = Vector3::new(
            f64::from(copy % 5) * SPACING,
            f64::from(copy / 5 % 5) * SPACING,
            f64::from(copy / 25) * SPACING,
        );
        let mut shifted = kras.clone();
        shifted.apply_rigid_transform(&RigidTransform::new(
            Matrix3::identity(),
            Quantity::new(offset, ANGSTROM),
        )?)?;
        editor.append_model(&shifted)?;
    }
    let model = editor.finish()?;
    let topology = model.shared_topology();

    let started = Instant::now();
    let features = InteractionFeatures::new(Arc::clone(&topology));
    let perceived = started.elapsed();
    println!(
        "informational only; target={}-{}; atoms={}; median of {PASSES} passes",
        std::env::consts::OS,
        std::env::consts::ARCH,
        model.atom_count(),
    );
    println!(
        "features perceived once in {:.1} ms",
        perceived.as_secs_f64() * 1e3
    );
    println!("| case | first atoms | second atoms | interactions | µs per search |");
    println!("|---|---:|---:|---:|---:|");

    // The inhibitor of a copy in the middle of the block, by atom order.
    let middle = (2 + 5 * 2 + 25) as usize * kras.atom_count();
    let ligand = AtomSelection::from_predicate(Arc::clone(&topology), |atom| {
        let index = atom.index().index();
        (middle..middle + kras.atom_count()).contains(&index)
            && atom
                .residue()
                .is_some_and(|residue| residue.name() == "A1L9E")
    });
    let all = AtomSelection::all(Arc::clone(&topology));
    let criteria = InteractionCriteria::default();
    let within = |first: &AtomSelection| {
        let mut found = None;
        let time = median(|| {
            found = Some(
                find_interactions_within(&features, &model, first, &all, &criteria)
                    .expect("valid search"),
            );
        });
        (time, found.expect("searched"))
    };
    let row = |case: &str, first: usize, second: usize, (time, found): (Duration, Interactions)| {
        println!(
            "| {case} | {first} | {second} | {} | {:.0} |",
            found.interactions().len(),
            time.as_secs_f64() * 1e6
        );
    };

    row(
        "inhibitor × every atom",
        ligand.len(),
        all.len(),
        within(&ligand),
    );
    let pocket = ligand.expand_to_residues();
    let marked = AtomSelection::from_predicate(Arc::clone(&topology), |atom| {
        let point = model.positions().values().value()[atom.index().index()];
        ligand
            .indices()
            .iter()
            .any(|other| (model.positions().values().value()[other.index()] - point).norm() <= 0.8)
    })
    .expand_to_residues()
    .union(&pocket)?;
    row(
        "inhibitor and pocket × every atom",
        marked.len(),
        all.len(),
        within(&marked),
    );
    let copy = AtomSelection::from_predicate(Arc::clone(&topology), |atom| {
        (middle..middle + kras.atom_count()).contains(&atom.index().index())
    });
    row(
        "one copy × every atom",
        copy.len(),
        all.len(),
        within(&copy),
    );

    // The inhibitor as its own structure against the block.
    let alone = model.subset(&ligand)?;
    let alone_features = InteractionFeatures::new(alone.shared_topology());
    let alone_atoms = AtomSelection::all(alone.shared_topology());
    // The block without that inhibitor, which the separate structure
    // replaces.
    let rest = all.difference(&ligand)?;
    let mut found = None;
    let time = median(|| {
        found = Some(
            find_interactions_between(
                &InteractionSide::new(&alone_features, &alone, &alone_atoms).expect("bound"),
                &InteractionSide::new(&features, &model, &rest).expect("bound"),
                &criteria,
            )
            .expect("valid search"),
        );
    });
    row(
        "inhibitor structure × block structure",
        alone_atoms.len(),
        rest.len(),
        (time, found.expect("searched")),
    );
    Ok(())
}

fn median(mut run: impl FnMut()) -> Duration {
    let mut times = (0..PASSES)
        .map(|_| {
            let started = Instant::now();
            run();
            black_box(started.elapsed())
        })
        .collect::<Vec<_>>();
    times.sort_unstable();
    times[PASSES / 2]
}
