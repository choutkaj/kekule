#![no_main]
use kekule::traj::frames::{read_frame_table_from, FrameTableWriter};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(table) = read_frame_table_from(data) else {
        return;
    };
    // An accepted table writes under its own header, and the written rows
    // read back to the same steps and cells.
    let header = table.header().clone();
    let mut writer = FrameTableWriter::new(Vec::new(), header.clone()).expect("a header writes");
    let mut row = writer.row();
    for (index, &step) in table.steps().iter().enumerate() {
        row.clear();
        row.set_step(step);
        for spec in header.columns() {
            let column = table.column(spec.key()).expect("a header column");
            if let Some(value) = column.value(index).expect("a row of the table") {
                row.set(spec.key(), value).expect("a read cell sets");
            }
        }
        writer.write_row(&row).expect("a read row writes");
    }
    let bytes = writer.into_inner().expect("a buffer flushes");
    let again = read_frame_table_from(bytes.as_slice()).expect("a written table reads");
    assert_eq!(again.header(), &header);
    assert_eq!(again.steps(), table.steps());
    assert_eq!(again.pending_bytes(), 0);
    for spec in header.columns() {
        assert_eq!(again.column(spec.key()), table.column(spec.key()));
    }
});
