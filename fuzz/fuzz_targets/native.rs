#![no_main]
use kekule::native::{self, NativeObject, NativeReadOptions, NativeWriteOptions};
use libfuzzer_sys::fuzz_target;

fn encode(object: &NativeObject) -> Vec<u8> {
    let options = NativeWriteOptions::default();
    match object {
        NativeObject::Molecule(value) => native::write(value.as_ref(), options),
        NativeObject::Topology(value) => native::write(value, options),
        NativeObject::Model(value) => native::write(value, options),
        NativeObject::Ensemble(value) => native::write(value, options),
        NativeObject::Trajectory(value) => native::write(value, options),
        NativeObject::Series(value) => native::write(value, options),
        _ => return Vec::new(),
    }
    .expect("a decoded object encodes")
}

fuzz_target!(|data: &[u8]| {
    let options = NativeReadOptions::default()
        .with_max_input_bytes(1 << 16)
        .expect("positive fuzz bounds")
        .with_max_atoms(4096)
        .expect("positive fuzz bounds")
        .with_max_bonds(8192)
        .expect("positive fuzz bounds")
        .with_max_items(256)
        .expect("positive fuzz bounds")
        .with_max_string_bytes(256)
        .expect("positive fuzz bounds");
    let mut input = data;
    let Ok(object) = native::read(&mut input, options) else {
        return;
    };
    // A decoded object is valid: it encodes, decodes again, and reaches a
    // fixed point, since equal values encode to equal bytes.
    let bytes = encode(&object);
    let again = native::read(bytes.as_slice(), NativeReadOptions::default())
        .expect("an encoded object decodes");
    assert!(encode(&again) == bytes, "re-encoding is not a fixed point");
});
