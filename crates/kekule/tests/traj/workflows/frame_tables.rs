//! Frame tables: exact round trips, crash and NUL tails, corrupt rows,
//! continuing at a checkpoint, following a growing file, and joining a
//! table into sets in memory and in files by step.

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use kekule::math::Point3;
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{
    ColumnRef, Geometry, Order, Positions, Sampling, StructureSet, StructureSetError,
    ThermodynamicState,
};
use kekule::traj::frames::{
    append_frame_table, create_frame_table, read_frame_table, read_frame_table_from, ColumnSpec,
    ColumnType, FrameTable, FrameTableErrorKind, FrameTableHeader, FrameTableReader,
    FrameTableWriter, KeepRows,
};
use kekule::traj::io::trr::{TrrLambdaPolicy, TrrWriteOptions};
use kekule::traj::io::{
    open_indexed_structure_set, write_structure_set, TrajectoryOpenOptions, TrajectoryWriteOptions,
};
use kekule::traj::TrajectoryFormat;
use kekule::units::{
    Quantity, BAR, DIMENSIONLESS, KELVIN, KILOCALORIE_PER_MOLE, KILOJOULE_PER_MOLE, NANOMETER,
    PICOSECOND,
};

use crate::support::{linear_carbon_topology, time_ordered, Frame};

/// A scratch directory removed on drop.
struct Directory(PathBuf);

impl Directory {
    fn new() -> Self {
        static SEQUENCE: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "kekule-frame-tables-{}-{}",
            std::process::id(),
            SEQUENCE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for Directory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

fn kj(value: f64) -> PropertyValue {
    PropertyValue::Real {
        value,
        unit: KILOJOULE_PER_MOLE,
    }
}

/// A GaMD-like header: time, stage, energies, a flag, a count, and the
/// sampling record of an NpT run with two bias terms.
fn header() -> FrameTableHeader {
    let state = ThermodynamicState::isothermal_isobaric(
        Quantity::new(300.0, KELVIN),
        Quantity::new(1.0, BAR),
    )
    .unwrap();
    let sampling = Sampling::new(state)
        .with_from_step(2_000)
        .with_method("GaMD dual boost\twith\\tabs")
        .with_potential_energy(key("potential_energy"))
        .unwrap()
        .with_bias([key("bias.potential"), key("bias.dihedral")])
        .unwrap();
    FrameTableHeader::new(vec![
        ColumnSpec::new(key("time"), ColumnType::Real(PICOSECOND)),
        ColumnSpec::new(key("stage"), ColumnType::Text),
        ColumnSpec::new(
            key("potential_energy"),
            ColumnType::Real(KILOJOULE_PER_MOLE),
        ),
        ColumnSpec::new(key("bias.potential"), ColumnType::Real(KILOJOULE_PER_MOLE)),
        ColumnSpec::new(key("bias.dihedral"), ColumnType::Real(KILOJOULE_PER_MOLE)),
        ColumnSpec::new(key("gamd.potential.k0"), ColumnType::Real(DIMENSIONLESS)),
        ColumnSpec::new(key("boosted"), ColumnType::Bool),
        ColumnSpec::new(key("replica"), ColumnType::Int),
    ])
    .unwrap()
    .with_time(key("time"))
    .unwrap()
    .with_sampling(sampling)
    .unwrap()
}

/// Writes rows of steps `steps`, 500 apart from 0, to `writer`.
fn write_rows<W: Write>(writer: &mut FrameTableWriter<W>, steps: impl IntoIterator<Item = u64>) {
    let mut row = writer.row();
    for index in steps {
        let step = 500 * index;
        row.clear();
        row.set_step(step);
        let x = index as f64;
        let set = |row: &mut kekule::traj::frames::FrameRow, name: &str, value| {
            row.set(&key(name), value).unwrap();
        };
        set(
            &mut row,
            "time",
            PropertyValue::Real {
                value: 0.002 * step as f64,
                unit: PICOSECOND,
            },
        );
        set(
            &mut row,
            "stage",
            PropertyValue::String(if index < 4 { "nvt" } else { "md" }.into()),
        );
        set(
            &mut row,
            "potential_energy",
            kj(-30_521.441_098_1 + x / 3.0),
        );
        set(&mut row, "bias.potential", kj(0.1 * x));
        if index % 3 != 1 {
            set(&mut row, "bias.dihedral", kj(1e-300 * x));
        }
        set(&mut row, "boosted", PropertyValue::Bool(index >= 4));
        set(&mut row, "replica", PropertyValue::Int(-(index as i64)));
        writer.write_row(&row).unwrap();
    }
}

fn real_values(table: &FrameTable, name: &str) -> Vec<Option<f64>> {
    match table.column(&key(name)).unwrap() {
        PropertyColumn::Real { values, .. } => values.clone(),
        other => panic!("{name} is not real: {other:?}"),
    }
}

#[test]
fn a_table_round_trips_its_header_and_every_bit_of_its_values() {
    let tricky = [
        -0.0,
        0.0,
        f64::MIN_POSITIVE,
        5e-324,
        -1e-300,
        1e300,
        -1.797_693_134_862_315_7e308,
        1.0 / 3.0,
        1e-4,
        9.999_999_999_999_999e-5,
        1e16,
        9_999_999_999_999_998.0,
        -30_521.441_098_1,
    ];
    let energy = key("energy");
    let label = key("label");
    let header = FrameTableHeader::new(vec![
        ColumnSpec::new(energy.clone(), ColumnType::Real(KILOJOULE_PER_MOLE)),
        ColumnSpec::new(label.clone(), ColumnType::Text),
    ])
    .unwrap();
    let texts = [
        Some("\\N"),
        Some(""),
        Some("tab\there"),
        Some("line\nbreak\r"),
        Some("back\\slash"),
        None,
        Some("ünïcode"),
    ];
    let mut writer = FrameTableWriter::new(Vec::new(), header.clone()).unwrap();
    let mut row = writer.row();
    for (index, value) in tricky.iter().chain(&[f64::NAN, f64::INFINITY]).enumerate() {
        row.clear();
        row.set_step(index as u64);
        row.set(&energy, kj(*value)).unwrap();
        if let Some(text) = texts.get(index).copied().flatten() {
            row.set(&label, PropertyValue::String(text.into())).unwrap();
        }
        writer.write_row(&row).unwrap();
    }
    let bytes = writer.into_inner().unwrap();
    let table = read_frame_table_from(bytes.as_slice()).unwrap();
    assert_eq!(table.header(), &header);
    let values = real_values(&table, "energy");
    for (read, written) in values.iter().zip(&tricky) {
        assert_eq!(read.unwrap().to_bits(), written.to_bits(), "{written:e}");
    }
    // Values that are not finite are written as missing.
    assert_eq!(&values[tricky.len()..], [None, None]);
    let PropertyColumn::String(labels) = table.column(&label).unwrap() else {
        panic!("text columns read as strings");
    };
    let expected = (0..values.len())
        .map(|index| texts.get(index).copied().flatten().map(str::to_owned))
        .collect::<Vec<_>>();
    assert_eq!(labels, &expected);
    let text = String::from_utf8(bytes).unwrap();
    assert!(text.contains("\n0\t-0\t\\\\N\n"), "{text}");
    assert!(text.contains("\t5e-324\t"), "{text}");
    assert!(text.contains("\n5\t1e300\t\\N\n"), "{text}");
}

#[test]
fn the_header_records_the_sampling_record_and_names_the_columns() {
    let mut writer = FrameTableWriter::new(Vec::new(), header()).unwrap();
    write_rows(&mut writer, 0..3);
    let bytes = writer.into_inner().unwrap();
    let text = String::from_utf8(bytes.clone()).unwrap();
    let lines = text.lines().collect::<Vec<_>>();
    assert_eq!(lines[0], "#kekule-frames\t1");
    assert_eq!(lines[1], "#column\tstep\tint");
    assert!(lines.contains(&"#column\tgamd.potential.k0\treal\t1"));
    assert!(lines.contains(&"#time\ttime"));
    assert!(lines.contains(&"#state\ttemperature\t300\tK"));
    assert!(lines.contains(&"#state\tpressure\t1\tbar"));
    assert!(lines.contains(&"#sampling\tfrom_step\t2000"));
    assert!(lines.contains(&"#sampling\tmethod\tGaMD dual boost\\twith\\\\tabs"));
    assert!(lines.contains(&"#bias\tbias.dihedral"));
    assert_eq!(
        lines[lines.len() - 4],
        "step\ttime\tstage\tpotential_energy\tbias.potential\tbias.dihedral\tgamd.potential.k0\tboosted\treplica"
    );
    let table = read_frame_table_from(bytes.as_slice()).unwrap();
    assert_eq!(table.header(), &header());
    assert_eq!(
        table.header().sampling().unwrap().method(),
        Some("GaMD dual boost\twith\\tabs")
    );
    assert_eq!(table.steps(), [0, 500, 1000]);
    // A missing real is an empty cell.
    assert_eq!(real_values(&table, "bias.dihedral")[1], None);

    // Comments in the header are ignored.
    let commented = text.replacen("#time", "## a note\n#time", 1);
    assert_eq!(read_frame_table_from(commented.as_bytes()).unwrap(), table);
}

#[test]
fn readers_leave_out_an_unfinished_line_and_nul_blocks() {
    let mut writer = FrameTableWriter::new(Vec::new(), header()).unwrap();
    write_rows(&mut writer, 0..3);
    let complete = writer.into_inner().unwrap();
    let table = read_frame_table_from(complete.as_slice()).unwrap();
    assert_eq!(table.pending_bytes(), 0);

    let mut crashed = complete.clone();
    crashed.extend_from_slice(b"1500\t3\tm");
    let read = read_frame_table_from(crashed.as_slice()).unwrap();
    assert_eq!(read.steps(), table.steps());
    assert_eq!(read.pending_bytes(), 8);

    let mut lost = complete.clone();
    lost.extend_from_slice(&[0; 4096]);
    let read = read_frame_table_from(lost.as_slice()).unwrap();
    assert_eq!(read, FrameTable::clone(&read));
    assert_eq!(read.steps(), table.steps());
    assert_eq!(read.pending_bytes(), 4096);

    // A header that is not complete is pending; a file of NUL blocks too.
    for bytes in [&complete[..40], &[0; 512][..], b""] {
        let error = read_frame_table_from(bytes).unwrap_err();
        assert_eq!(error.kind(), FrameTableErrorKind::Pending);
    }
    let error = read_frame_table_from(b"step\tbias".as_slice()).unwrap_err();
    assert_eq!(error.kind(), FrameTableErrorKind::Corrupt);
}

#[test]
fn a_complete_row_that_does_not_parse_is_corrupt() {
    let mut writer = FrameTableWriter::new(Vec::new(), header()).unwrap();
    write_rows(&mut writer, 0..6);
    let text = String::from_utf8(writer.into_inner().unwrap()).unwrap();
    let names = text
        .lines()
        .position(|line| line.starts_with("step"))
        .unwrap()
        + 1;
    for (damage, line) in [
        (text.replacen("\n500\t1\t", "\n500\t1x\t", 1), names + 2),
        (text.replacen("\n500\t", "\n-500\t", 1), names + 2),
        (text.replacen("\n1000\t", "\n400\t", 1), names + 3),
        (text.replacen("\tnvt\t", "\tnvt\textra\t", 1), names + 1),
        (text.replacen("\ttrue\t", "\tyes\t", 1), names + 5),
        (text.replacen("\tnvt\t", "\tn\\xvt\t", 1), names + 1),
    ] {
        let error = read_frame_table_from(damage.as_bytes()).unwrap_err();
        assert_eq!(error.kind(), FrameTableErrorKind::Corrupt, "{error}");
        assert_eq!(error.line(), Some(line), "{error}");
    }
    let names_line = text.replacen("step\ttime", "step\tthyme", 1);
    let error = read_frame_table_from(names_line.as_bytes()).unwrap_err();
    assert_eq!(error.kind(), FrameTableErrorKind::Corrupt);
    assert_eq!(error.line(), Some(names));
}

#[test]
fn other_versions_and_unknown_records_are_unsupported() {
    let text = String::from_utf8(
        FrameTableWriter::new(Vec::new(), header())
            .unwrap()
            .into_inner()
            .unwrap(),
    )
    .unwrap();
    let newer = text.replacen("#kekule-frames\t1", "#kekule-frames\t2", 1);
    let error = read_frame_table_from(newer.as_bytes()).unwrap_err();
    assert_eq!(error.kind(), FrameTableErrorKind::Unsupported);
    assert_eq!(
        error.to_string(),
        "The frame table has version \"2\". Only version 1 is read."
    );
    let record = text.replacen("#time", "#sampling\tvolume\tvolume\n#time", 1);
    let error = read_frame_table_from(record.as_bytes()).unwrap_err();
    assert_eq!(error.kind(), FrameTableErrorKind::Unsupported);
    assert!(error.line().is_some());
}

#[test]
fn writers_and_headers_refuse_what_the_format_cannot_hold() {
    let error =
        FrameTableHeader::new(vec![ColumnSpec::new(key("step"), ColumnType::Int)]).unwrap_err();
    assert_eq!(error.kind(), FrameTableErrorKind::Invalid);
    let twice = || ColumnSpec::new(key("x"), ColumnType::Int);
    let error = FrameTableHeader::new(vec![twice(), twice()]).unwrap_err();
    assert_eq!(error.column(), Some(&key("x")));
    let unnamed = kekule::units::Unit::new(NANOMETER.dimension(), 2e-9, None).unwrap();
    let error = FrameTableHeader::new(vec![ColumnSpec::new(key("x"), ColumnType::Real(unnamed))])
        .unwrap_err();
    assert_eq!(error.column(), Some(&key("x")));
    let energies = FrameTableHeader::new(vec![
        ColumnSpec::new(key("e"), ColumnType::Real(KILOJOULE_PER_MOLE)),
        ColumnSpec::new(key("d"), ColumnType::Real(NANOMETER)),
    ])
    .unwrap();
    assert!(energies.clone().with_time(key("e")).is_err());
    assert!(energies.clone().with_time(key("t")).is_err());
    let state = ThermodynamicState::canonical(Quantity::new(300.0, KELVIN)).unwrap();
    let biased = |name: &str| Sampling::new(state).with_bias([key(name)]).unwrap();
    assert_eq!(
        energies
            .clone()
            .with_sampling(biased("d"))
            .unwrap_err()
            .column(),
        Some(&key("d"))
    );
    assert!(energies.clone().with_sampling(biased("q")).is_err());

    let mut writer = FrameTableWriter::new(Vec::new(), energies.clone()).unwrap();
    let mut row = writer.row();
    assert_eq!(
        row.set(&key("d"), kj(1.0)).unwrap_err().kind(),
        FrameTableErrorKind::Invalid
    );
    assert!(row.set(&key("e"), PropertyValue::Int(1)).is_err());
    assert!(row.set(&key("q"), kj(1.0)).is_err());
    // A value converts to its column's unit.
    row.set(
        &key("e"),
        PropertyValue::Real {
            value: 1.0,
            unit: KILOCALORIE_PER_MOLE,
        },
    )
    .unwrap();
    assert_eq!(row.value(&key("e")), Some(&kj(4.184)));
    assert!(writer.write_row(&row).is_err(), "a row needs a step");
    row.set_step(10);
    writer.write_row(&row).unwrap();
    let error = writer.write_row(&row).unwrap_err();
    assert_eq!(error.to_string(), "Step 10 does not follow step 10.");
    let other = FrameTableWriter::new(Vec::new(), header()).unwrap().row();
    assert!(writer.write_row(&other).is_err());
    assert_eq!(writer.last_step(), Some(10));
}

#[test]
fn a_table_continued_at_a_checkpoint_equals_one_written_straight_through() {
    let directory = Directory::new();
    let straight = directory.path("straight.tsv");
    let mut writer = create_frame_table(&straight, header()).unwrap();
    write_rows(&mut writer, 0..10);
    drop(writer);

    let continued = directory.path("continued.tsv");
    let mut writer = create_frame_table(&continued, header()).unwrap();
    write_rows(&mut writer, 0..7);
    drop(writer);
    // The run died while writing a row.
    let mut file = fs::OpenOptions::new()
        .append(true)
        .open(&continued)
        .unwrap();
    file.write_all(b"3500\t7\tm").unwrap();
    file.write_all(&[0; 100]).unwrap();
    drop(file);

    let (mut writer, report) =
        append_frame_table(&continued, &header(), KeepRows::ThroughStep(2_000)).unwrap();
    assert_eq!(report.kept_rows(), 5);
    assert_eq!(report.last_step(), Some(2_000));
    assert_eq!(report.dropped_rows(), 2);
    assert!(report.dropped_bytes() > 108);
    assert_eq!(writer.last_step(), Some(2_000));
    write_rows(&mut writer, 5..10);
    drop(writer);
    assert_eq!(fs::read(&continued).unwrap(), fs::read(&straight).unwrap());

    // Keeping every row cuts only the pending tail.
    let mut file = fs::OpenOptions::new()
        .append(true)
        .open(&continued)
        .unwrap();
    file.write_all(b"5000\t").unwrap();
    drop(file);
    let (writer, report) = append_frame_table(&continued, &header(), KeepRows::All).unwrap();
    assert_eq!((report.kept_rows(), report.dropped_bytes()), (10, 5));
    drop(writer);
    assert_eq!(fs::read(&continued).unwrap(), fs::read(&straight).unwrap());

    // A checkpoint before the first row keeps the header alone.
    let (mut writer, report) =
        append_frame_table(&continued, &header(), KeepRows::ThroughStep(0)).unwrap();
    assert_eq!(report.kept_rows(), 1);
    write_rows(&mut writer, 1..10);
    drop(writer);
    assert_eq!(fs::read(&continued).unwrap(), fs::read(&straight).unwrap());
}

#[test]
fn a_header_that_differs_leaves_the_file_unchanged() {
    let directory = Directory::new();
    let path = directory.path("frames.tsv");
    let mut writer = create_frame_table(&path, header()).unwrap();
    write_rows(&mut writer, 0..3);
    drop(writer);
    let mut file = fs::OpenOptions::new().append(true).open(&path).unwrap();
    file.write_all(b"1500\t").unwrap();
    drop(file);
    let before = fs::read(&path).unwrap();
    let state = ThermodynamicState::canonical(Quantity::new(310.0, KELVIN)).unwrap();
    let other = FrameTableHeader::new(header().columns().to_vec())
        .unwrap()
        .with_sampling(Sampling::new(state))
        .unwrap();
    let error = append_frame_table(&path, &other, KeepRows::ThroughStep(500)).unwrap_err();
    assert_eq!(error.kind(), FrameTableErrorKind::Mismatch);
    assert_eq!(fs::read(&path).unwrap(), before);

    // A corrupt row leaves it unchanged too.
    let corrupt = String::from_utf8(before.clone())
        .unwrap()
        .replacen("\n500\t", "\n500x\t", 1);
    fs::write(&path, &corrupt).unwrap();
    let error = append_frame_table(&path, &header(), KeepRows::All).unwrap_err();
    assert_eq!(error.kind(), FrameTableErrorKind::Corrupt);
    assert_eq!(fs::read(&path).unwrap(), corrupt.as_bytes());
}

#[test]
fn a_table_without_a_complete_header_is_started_again() {
    let directory = Directory::new();
    let path = directory.path("frames.tsv");
    let straight = FrameTableWriter::new(Vec::new(), header())
        .unwrap()
        .into_inner()
        .unwrap();
    fs::write(&path, &straight[..30]).unwrap();
    let (writer, report) = append_frame_table(&path, &header(), KeepRows::ThroughStep(0)).unwrap();
    assert_eq!(report.kept_rows(), 0);
    assert_eq!(report.last_step(), None);
    assert_eq!(report.dropped_bytes(), 30);
    drop(writer);
    assert_eq!(fs::read(&path).unwrap(), straight);
}

#[test]
fn a_growing_table_is_read_row_by_row_until_it_is_cut() {
    let directory = Directory::new();
    let path = directory.path("frames.tsv");
    let mut writer = create_frame_table(&path, header()).unwrap();
    let mut reader = FrameTableReader::open(&path).unwrap();
    assert!(reader.table().is_empty());
    write_rows(&mut writer, 0..2);
    assert_eq!(reader.refresh().unwrap(), 2);
    assert_eq!(reader.refresh().unwrap(), 0);
    write_rows(&mut writer, 2..5);
    writer.get_ref().get_ref().write_all(b"2500\t5").unwrap();
    assert_eq!(reader.refresh().unwrap(), 3);
    assert_eq!(reader.table().pending_bytes(), 6);
    assert_eq!(reader.table(), &read_frame_table(&path).unwrap());
    drop(writer);

    // The run continues from its checkpoint at step 1000.
    let (mut writer, _) =
        append_frame_table(&path, &header(), KeepRows::ThroughStep(1_000)).unwrap();
    write_rows(&mut writer, 3..4);
    let before = reader.table().clone();
    let error = reader.refresh().unwrap_err();
    assert_eq!(error.kind(), FrameTableErrorKind::Changed);
    assert_eq!(reader.table(), &before);
    let reopened = FrameTableReader::open(&path).unwrap();
    assert_eq!(reopened.table().steps(), [0, 500, 1000, 1500]);
}

/// A table of `rows` rows 500 steps apart, with a bias in every row.
fn bias_table(rows: u64) -> FrameTable {
    let mut writer = FrameTableWriter::new(Vec::new(), header()).unwrap();
    write_rows(&mut writer, 0..rows);
    read_frame_table_from(writer.into_inner().unwrap().as_slice()).unwrap()
}

/// Geometries of a twelve-atom chain at steps `steps`.
fn frames(steps: &[u64]) -> StructureSet {
    let topology = linear_carbon_topology(12);
    time_ordered(
        Arc::clone(&topology),
        steps.iter().map(|&step| {
            let points = (0..12)
                .map(|atom| Point3::new(0.15 * atom as f64, 0.01 * step as f64 / 500.0, 0.0))
                .collect();
            let mut frame =
                Frame::new(Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap());
            frame.set_step(Some(step));
            frame
                .set_time(Some(Quantity::new(0.002 * step as f64, PICOSECOND)))
                .unwrap();
            frame
        }),
    )
}

#[test]
fn joining_a_table_fills_columns_by_step_and_installs_its_sampling_record() {
    let mut set = frames(&[0, 1_000, 2_000, 2_250]);
    let table = bias_table(5);
    let report = set.join_frames(&table).unwrap();
    assert_eq!((report.matched(), report.unmatched()), (3, 1));
    assert_eq!(
        report.columns(),
        [
            key("potential_energy"),
            key("bias.potential"),
            key("bias.dihedral"),
            key("gamd.potential.k0"),
            key("replica"),
        ]
    );
    // The time, text, and boolean columns stay in the table.
    for skipped in ["time", "stage", "boosted", "step"] {
        assert!(set.columns().get(&key(skipped)).is_none(), "{skipped}");
    }
    let bias = set
        .column_values(&ColumnRef::Key(key("bias.potential")))
        .unwrap();
    assert_eq!(bias.value(), &[Some(0.0), Some(0.2), Some(0.4), None]);
    assert_eq!(set.sampling(), table.header().sampling());
    // The k0 column is all missing in the table, so the set has none.
    assert!(set.columns().get(&key("gamd.potential.k0")).is_none());

    // A grown table fills the cells its new rows hold; joining is idempotent.
    let mut grown = set.clone();
    let mut longer = frames(&[0, 1_000, 2_000, 2_250, 2_500]);
    longer.join_frames(&bias_table(5)).unwrap();
    grown.join_frames(&bias_table(4)).unwrap();
    longer.join_frames(&bias_table(7)).unwrap();
    let again = longer.clone();
    longer.join_frames(&bias_table(7)).unwrap();
    assert_eq!(longer.columns(), again.columns());
    assert_eq!(
        longer
            .column_values(&ColumnRef::Key(key("bias.potential")))
            .unwrap()
            .value(),
        &[Some(0.0), Some(0.2), Some(0.4), None, Some(0.5)]
    );

    // A different sampling record is refused and changes nothing.
    let state = ThermodynamicState::canonical(Quantity::new(310.0, KELVIN)).unwrap();
    let other = FrameTableHeader::new(header().columns().to_vec())
        .unwrap()
        .with_sampling(Sampling::new(state))
        .unwrap();
    let mut writer = FrameTableWriter::new(Vec::new(), other).unwrap();
    write_rows(&mut writer, 0..2);
    let other = read_frame_table_from(writer.into_inner().unwrap().as_slice()).unwrap();
    let before = set.clone();
    assert_eq!(
        set.join_frames(&other).unwrap_err(),
        StructureSetError::SamplingMismatch
    );
    assert_eq!(set.columns(), before.columns());
    assert_eq!(set.sampling(), before.sampling());
}

#[test]
fn a_set_without_order_steps_joins_by_its_step_column() {
    let topology = linear_carbon_topology(12);
    let geometry = || Geometry::new(Positions::zeros(topology.atom_count()));
    let mut set =
        StructureSet::from_geometries(Arc::clone(&topology), [geometry(), geometry()]).unwrap();
    assert_eq!(
        set.join_frames(&bias_table(3)).unwrap_err(),
        StructureSetError::Missing(kekule::structure::Requirement::Steps)
    );
    set.columns_mut()
        .insert(
            key("step"),
            PropertyColumn::Int(vec![Some(1_000), Some(-1)]),
        )
        .unwrap();
    let report = set.join_frames(&bias_table(3)).unwrap();
    assert_eq!((report.matched(), report.unmatched()), (1, 1));
    set.set_order(Some(Order::Step(vec![0, 500]))).unwrap();
    set.join_frames(&bias_table(3)).unwrap();
    assert_eq!(
        set.column_values(&ColumnRef::Key(key("bias.potential")))
            .unwrap()
            .value(),
        &[Some(0.0), Some(0.1)]
    );
}

#[test]
fn a_set_read_from_a_file_joins_the_table_and_keeps_it_on_refresh() {
    let directory = Directory::new();
    let path = directory.path("run.trr");
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Trr).with_trr_options(
        TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
    );
    write_structure_set(&path, &frames(&[0, 1_000, 2_000]), options).unwrap();
    let topology = linear_carbon_topology(12);
    let growing = TrajectoryOpenOptions::default().with_growing(true);
    let mut set = open_indexed_structure_set(&path, topology, growing).unwrap();
    let report = set.join_frames(&bias_table(4)).unwrap();
    assert_eq!(report.matched(), 2);
    assert_eq!(
        set.column_values(&ColumnRef::Key(key("potential_energy")))
            .unwrap()
            .value()
            .iter()
            .map(Option::is_some)
            .collect::<Vec<_>>(),
        [true, true, false]
    );
    let sampling = set.sampling().cloned();
    assert!(sampling.is_some());
    let loaded = StructureSet::from_reader(&set, [2, 0]).unwrap();
    assert_eq!(loaded.sampling(), sampling.as_ref());
    set.join_frames(&bias_table(5)).unwrap();
    assert_eq!(
        set.column_values(&ColumnRef::Key(key("bias.potential")))
            .unwrap()
            .value(),
        &[Some(0.0), Some(0.2), Some(0.4)]
    );
    assert_eq!(set.refresh().unwrap(), 0);
    assert_eq!(set.sampling(), sampling.as_ref());
}
