//! Optional comparison against independently generated periodic references.
//!
//! Run with TOPOLOGY.txt INPUT.xtc REFERENCES_DIR. Generate references from
//! externally supplied data using benchmarks/reference/trajectory/export_periodic.py.
//! This is a scientific development check, not a routine CI or release gate.

use kekule::{
    core::{Atom, BondOrder, Element, MoleculeEditor},
    topology::{AtomSelection, Topology},
    units::NANOMETER,
};
use kekule_traj::{
    io::read_trajectory,
    periodic::{MoleculeImager, TrajectoryUnwrapper},
    FrameBuffer, Trajectory,
};
use std::{error::Error, fs, path::Path, sync::Arc};

fn topology(path: &Path) -> Result<Arc<Topology>, Box<dyn Error>> {
    let input = fs::read_to_string(path)?;
    let mut words = input.split_whitespace();
    let atoms: usize = words.next().ok_or("missing atom count")?.parse()?;
    let bonds: usize = words.next().ok_or("missing bond count")?.parse()?;
    let mut molecule = MoleculeEditor::new();
    let mut ids = Vec::with_capacity(atoms);
    for _ in 0..atoms {
        let element = Element::from_symbol(words.next().ok_or("missing element")?)
            .ok_or("invalid element")?;
        ids.push(molecule.add_atom(Atom::new(element))?);
    }
    for _ in 0..bonds {
        let a: usize = words.next().ok_or("missing bond endpoint")?.parse()?;
        let b: usize = words.next().ok_or("missing bond endpoint")?.parse()?;
        // Only externally supplied connectivity matters to these coordinate operations.
        molecule.add_bond(
            *ids.get(a).ok_or("invalid atom")?,
            *ids.get(b).ok_or("invalid atom")?,
            BondOrder::Single,
        )?;
    }
    if words.next().is_some() {
        return Err("unexpected topology data".into());
    }
    Ok(Arc::new(Topology::from_molecule(&molecule.finish()?)?))
}

fn compare(actual: &Trajectory, path: &Path) -> Result<f64, Box<dyn Error>> {
    let reference = fs::read_to_string(path)?
        .split_whitespace()
        .map(str::parse::<f64>)
        .collect::<Result<Vec<_>, _>>()?;
    if reference.len() != actual.len() * actual.topology().atom_count() * 3 {
        return Err("reference coordinate dimensions differ".into());
    }
    let coordinates = actual
        .frames()
        .flat_map(|frame| frame.positions().values().value().to_vec())
        .flat_map(|point| [point.x, point.y, point.z]);
    let mut maximum = 0.0_f64;
    for (a, b) in coordinates.zip(reference) {
        if !a.is_finite() || !b.is_finite() {
            return Err("non-finite reference comparison".into());
        }
        maximum = maximum.max((a - b).abs());
    }
    // Predeclared tolerance accommodates independent f32 arithmetic, and is
    // smaller than the source XTC's usual 0.001 nm coordinate resolution.
    if maximum > 1.0e-4 {
        return Err(format!(
            "{}: maximum deviation {maximum} nm exceeds 0.0001 nm",
            path.display()
        )
        .into());
    }
    Ok(maximum)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = std::env::args_os().skip(1).collect::<Vec<_>>();
    if args.len() != 3 {
        return Err(
            "usage: trajectory_periodic_reference TOPOLOGY.txt INPUT.xtc REFERENCES_DIR".into(),
        );
    }
    let topology = topology(Path::new(&args[0]))?;
    let source = read_trajectory(Path::new(&args[1]), topology.clone())?;
    let anchors = AtomSelection::all(&topology);
    let whole = source.make_molecules_whole()?;
    let imaged = source.image_molecules(&anchors)?;
    let unwrapped = source.unwrap()?;
    let imager = MoleculeImager::new(topology.clone());
    let mut unwrapper = TrajectoryUnwrapper::new(topology.clone());
    let mut buffer = FrameBuffer::new(topology.clone());
    for (index, frame) in source.frames().enumerate() {
        let expected = whole.frame(index).unwrap();
        assert_eq!(imager.make_whole(index, frame)?, expected.to_frame());
        assert_eq!(
            imager.image(index, frame, &anchors)?,
            imaged.frame(index).unwrap().to_frame()
        );
        buffer.copy_from(frame)?;
        unwrapper.unwrap_in_place(index, &mut buffer)?;
        assert_eq!(
            buffer.frame_view().to_frame(),
            unwrapped.frame(index).unwrap().to_frame()
        );
        // Every non-position field survives each operation unchanged.
        for output in [
            expected,
            imaged.frame(index).unwrap(),
            unwrapped.frame(index).unwrap(),
        ] {
            assert_eq!(output.properties(), frame.properties());
            assert_eq!(output.cell(), frame.cell());
            assert_eq!(output.velocities(), frame.velocities());
            assert_eq!(output.forces(), frame.forces());
            assert_eq!(output.time(), frame.time());
            assert_eq!(output.step(), frame.step());
            assert_eq!(output.positions().values().unit(), NANOMETER);
        }
    }
    println!(
        "{{\"frames\":{},\"atoms\":{},\"maximum_error_nm\":{{",
        source.len(),
        topology.atom_count()
    );
    for (index, (label, trajectory)) in [
        ("raw", &source),
        ("whole", &whole),
        ("image", &imaged),
        ("unwrap", &unwrapped),
    ]
    .into_iter()
    .enumerate()
    {
        let error = compare(
            trajectory,
            &Path::new(&args[2]).join(format!("{label}.txt")),
        )?;
        println!("{}\"{label}\":{error}", if index == 0 { "" } else { "," });
    }
    println!("}},\"streaming_parity\":true,\"metadata_preserved\":true}}");
    Ok(())
}
