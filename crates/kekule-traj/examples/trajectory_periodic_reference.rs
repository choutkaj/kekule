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
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{error::Error, fs, io::Read, path::Path, sync::Arc};

fn sha256(path: &Path) -> Result<String, Box<dyn Error>> {
    let mut file = fs::File::open(path)?;
    let mut hash = Sha256::new();
    let mut buffer = [0; 65536];
    loop {
        let n = file.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        hash.update(&buffer[..n]);
    }
    Ok(hash
        .finalize()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect())
}
fn provenance(
    topology: &Path,
    trajectory: &Path,
    directory: &Path,
) -> Result<Value, Box<dyn Error>> {
    let metadata: Value =
        serde_json::from_reader(fs::File::open(directory.join("provenance.json"))?)?;
    if metadata["schema"] != 2
        || metadata["coordinate_unit"] != "nm"
        || metadata["tolerance_nm"] != json!(0.0001)
        || metadata["references"]["mdtraj"] != "1.11.1.post1"
        || metadata["references"]["MDAnalysis"] != "2.9.0"
    {
        return Err("incompatible reference provenance".into());
    }
    for (path, expected) in [
        (topology, &metadata["artifacts"]["topology.txt"]),
        (trajectory, &metadata["inputs"]["trajectory"]["sha256"]),
    ] {
        if expected.as_str() != Some(sha256(path)?.as_str()) {
            return Err("reference input checksum differs".into());
        }
    }
    for name in ["raw.txt", "whole.txt", "image.txt", "unwrap.txt"] {
        if metadata["artifacts"][name].as_str() != Some(sha256(&directory.join(name))?.as_str()) {
            return Err("reference artifact checksum differs".into());
        }
    }
    Ok(metadata)
}

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
    let provenance = provenance(
        Path::new(&args[0]),
        Path::new(&args[1]),
        Path::new(&args[2]),
    )?;
    let topology = topology(Path::new(&args[0]))?;
    let source = read_trajectory(Path::new(&args[1]), topology.clone())?;
    if source.is_empty()
        || topology.atom_count() == 0
        || provenance["frames"] != json!(source.len())
        || provenance["atoms"] != json!(topology.atom_count())
    {
        return Err("reference dimensions differ or are empty".into());
    }
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
    let mut errors = serde_json::Map::new();
    for (label, trajectory) in [
        ("raw", &source),
        ("whole", &whole),
        ("image", &imaged),
        ("unwrap", &unwrapped),
    ]
    .into_iter()
    {
        let error = compare(
            trajectory,
            &Path::new(&args[2]).join(format!("{label}.txt")),
        )?;
        errors.insert(label.into(), json!(error));
    }
    println!(
        "{}",
        json!({"frames":source.len(),"atoms":topology.atom_count(),"maximum_error_nm":errors,"streaming_parity":true,"metadata_preserved":true,"provenance":provenance})
    );
    Ok(())
}

#[test]
fn provenance_rejects_changed_sources_artifacts_and_contracts() {
    let directory = std::env::temp_dir().join(format!(
        "kekule-trajectory-provenance-{}",
        std::process::id()
    ));
    fs::create_dir(&directory).unwrap();
    let names = [
        "topology.txt",
        "raw.txt",
        "whole.txt",
        "image.txt",
        "unwrap.txt",
        "input.xtc",
    ];
    for name in names {
        fs::write(directory.join(name), b"regression bytes").unwrap();
    }
    let digest = sha256(&directory.join("input.xtc")).unwrap();
    let mut artifacts = serde_json::Map::new();
    for name in &names[..5] {
        artifacts.insert((*name).into(), json!(digest));
    }
    let baseline = json!({"schema":2,"coordinate_unit":"nm","tolerance_nm":0.0001,"references":{"mdtraj":"1.11.1.post1","MDAnalysis":"2.9.0"},"inputs":{"trajectory":{"sha256":digest}},"artifacts":artifacts});
    let write = |value: &Value| {
        fs::write(
            directory.join("provenance.json"),
            serde_json::to_vec(value).unwrap(),
        )
        .unwrap()
    };
    let check = || {
        provenance(
            &directory.join("topology.txt"),
            &directory.join("input.xtc"),
            &directory,
        )
    };
    write(&baseline);
    assert!(check().is_ok());
    for name in names {
        fs::write(directory.join(name), b"changed").unwrap();
        assert!(check().is_err(), "{name}");
        fs::write(directory.join(name), b"regression bytes").unwrap();
    }
    for (key, value) in [
        ("schema", json!(1)),
        ("coordinate_unit", json!("angstrom")),
        ("tolerance_nm", json!(0.1)),
    ] {
        let mut modified = baseline.clone();
        modified[key] = value;
        write(&modified);
        assert!(check().is_err());
    }
    for name in names {
        fs::remove_file(directory.join(name)).unwrap();
    }
    fs::remove_file(directory.join("provenance.json")).unwrap();
    fs::remove_dir(directory).unwrap();
}
