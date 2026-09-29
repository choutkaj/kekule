#![no_main]

use kekule::{perception::resonance::*, smiles};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    // External one-record seed files include a trailing line ending.
    let Ok(molecules) = smiles::to_molecules(input.trim()) else {
        return;
    };
    let mask = data.iter().fold(0u8, |sum, b| sum.wrapping_add(*b)) & 31;
    for mut mol in molecules {
        if mol.atom_count() > 32 || mol.bond_count() > 48 || mol.perceive().is_err() {
            continue;
        }
        let source = mol.clone();
        perceive_resonance(&mut mol).expect("default perception supplies conjugation");
        let before = mol.perception().clone();
        let options = ResonanceOptions {
            flags: ResonanceFlags::from_bits(mask).unwrap(),
            max_structures: 16,
            max_total_work: 10_000,
        };
        if let Ok(forms) = enumerate_resonance(&mol, options) {
            for i in 0..forms.contributors().len() {
                let form = forms
                    .to_molecule(i)
                    .expect("validated assignments materialize");
                assert_eq!(
                    form.atom_ids().collect::<Vec<_>>(),
                    source.atom_ids().collect::<Vec<_>>()
                );
                assert_eq!(
                    form.bond_ids().collect::<Vec<_>>(),
                    source.bond_ids().collect::<Vec<_>>()
                );
                assert!(!form.perception().has_conjugation());
            }
        }
        assert_eq!(mol, source);
        assert_eq!(mol.perception(), &before);
    }
});
