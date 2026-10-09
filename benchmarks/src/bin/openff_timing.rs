//! Repeatable warm CPU timings on externally supplied molecular inputs.
use kekule_openff::{ForceField, NaglModel};
use serde_json::{json, Value};
use std::{
    hint::black_box,
    io::{self, BufRead},
    time::Instant,
};

fn measure(
    mut operation: impl FnMut() -> Result<(), Box<dyn std::error::Error>>,
    count: usize,
) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    operation()?;
    (0..count)
        .map(|_| {
            let start = Instant::now();
            operation()?;
            Ok(start.elapsed().as_secs_f64() * 1000.0)
        })
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let model = NaglModel::load(
        std::env::args()
            .nth(1)
            .ok_or("expected model bundle path")?,
    )?;
    let model_load_ms = start.elapsed().as_secs_f64() * 1000.0;
    let ff = ForceField::rosemary()?;
    println!("{}", json!({"model_load_ms": model_load_ms}));
    for line in io::stdin().lock().lines() {
        let record: Value = serde_json::from_str(&line?)?;
        let mut molecules =
            kekule::smiles::to_molecules(record["smiles"].as_str().ok_or("missing smiles")?)?;
        if molecules.len() != 1 {
            return Err("expected one connected molecule".into());
        }
        let m = molecules.remove(0);
        let features = measure(
            || {
                black_box(kekule_openff::diagnostics::atom_features(&model, &m)?);
                Ok(())
            },
            5,
        )?;
        let infer = measure(
            || {
                black_box(kekule_openff::diagnostics::infer_charges(&model, &m)?);
                Ok(())
            },
            5,
        )?;
        let assign = measure(
            || {
                black_box(model.assign_charges(&m)?);
                Ok(())
            },
            5,
        )?;
        let full = measure(
            || {
                black_box(ff.parameterize_molecule(m.clone(), &model)?);
                Ok(())
            },
            3,
        )?;
        println!(
            "{}",
            json!({"id":record["id"],"atoms":m.atom_count(),"features_ms":features,
            "inference_including_features_ms":infer,"assign_charges_ms":assign,"full_parameterization_ms":full})
        );
    }
    Ok(())
}
