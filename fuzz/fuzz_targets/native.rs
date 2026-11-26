#![no_main]
use kekule::native::{self, NativeObject, NativeReadOptions, NativeWriteOptions};
use libfuzzer_sys::fuzz_target;

fn encode(object: &NativeObject) -> Option<Vec<u8>> {
    let options = NativeWriteOptions::default();
    let bytes = match object {
        NativeObject::Molecule(value) => native::write(value.as_ref(), options),
        NativeObject::Topology(value) => native::write(value, options),
        NativeObject::Structure(value) => native::write(value, options),
        NativeObject::StructureInterpretation(value) => native::write(value.as_ref(), options),
        NativeObject::StructureSet(value) => native::write(value, options),
        NativeObject::StructureSetInterpretation(value) => native::write(value.as_ref(), options),
        NativeObject::StructureSetDescription(value) => native::write(value.as_ref(), options),
        NativeObject::Weights(value) => native::write(value.as_ref(), options),
        NativeObject::Clustering(value) => native::write(value.as_ref(), options),
        _ => return None,
    };
    Some(bytes.expect("a decoded object encodes"))
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
    // A decoded object of any readable version is valid: it encodes in the
    // current version, decodes again, and reaches a fixed point, since equal
    // values encode to equal bytes.
    let Some(bytes) = encode(&object) else {
        return;
    };
    let again = native::read(bytes.as_slice(), NativeReadOptions::default())
        .expect("an encoded object decodes");
    assert!(
        encode(&again).as_ref() == Some(&bytes),
        "re-encoding is not a fixed point"
    );
});
