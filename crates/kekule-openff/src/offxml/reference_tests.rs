use super::*;
use serde_json::{json, Value};

fn fixture() -> Value {
    serde_json::from_reader(flate2::read::GzDecoder::new(
        &include_bytes!("../../tests/fixtures/offxml-generalization.json.gz")[..],
    ))
    .unwrap()
}

fn snapshot(ff: &ForceField) -> Value {
    let torsion = |rules: &[Rule<TorsionParameter>]| {
        rules.iter().map(|r| {
            let p = &r.parameter;
            json!({"id":p.source.id,"smirks":p.source.smirks,
                "terms":p.terms.iter().map(|t|json!({"k":t.k.value(),"phase":t.phase.value(),"periodicity":t.periodicity,"idivf":t.idivf})).collect::<Vec<_>>()})
        }).collect::<Vec<_>>()
    };
    let s = &ff.settings;
    json!({
        "parameters": {
            "Bonds":ff.bonds.iter().map(|r| {let p=&r.parameter;json!({"id":p.source.id,"smirks":p.source.smirks,"length":p.length.value(),"k":p.k.value()})}).collect::<Vec<_>>(),
            "Angles":ff.angles.iter().map(|r| {let p=&r.parameter;json!({"id":p.source.id,"smirks":p.source.smirks,"angle":p.angle.value(),"k":p.k.value()})}).collect::<Vec<_>>(),
            "ProperTorsions":torsion(&ff.propers),
            "ImproperTorsions":torsion(&ff.impropers),
            "Constraints":ff.constraints.iter().map(|r| {let p=&r.parameter;json!({"id":p.0.id,"smirks":p.0.smirks,"distance":p.1.as_ref().map(Quantity::value)})}).collect::<Vec<_>>(),
            "vdW":ff.vdw.iter().map(|r| {let p=&r.parameter;json!({"id":p.source.id,"smirks":p.source.smirks,"sigma":p.sigma.value(),"epsilon":p.epsilon.value()})}).collect::<Vec<_>>(),
            "LibraryCharges":ff.library.iter().map(|r| {let p=&r.parameter;json!({"id":p.source.id,"smirks":p.source.smirks,"charges":p.charges})}).collect::<Vec<_>>()
        },
        "settings":{
            "vdw_cutoff":s.vdw_cutoff.value(), "vdw_switch_width":s.vdw_switch_width.value(),
            "electrostatics_cutoff":s.electrostatics_cutoff.value(),"electrostatics_switch_width":s.electrostatics_switch_width.value(),
            "vdw_scales":s.vdw_scales,"electrostatics_scales":s.electrostatics_scales,
            "vdw_periodic_method":s.vdw_periodic_method,"vdw_nonperiodic_method":s.vdw_nonperiodic_method,
            "electrostatics_periodic_method":s.electrostatics_periodic_method,"electrostatics_nonperiodic_method":s.electrostatics_nonperiodic_method
        }
    })
}

fn compare(actual: &Value, expected: &Value, context: &str) {
    match (actual, expected) {
        (Value::Number(a), Value::Number(b)) => {
            let (a, b) = (a.as_f64().unwrap(), b.as_f64().unwrap());
            assert!(
                (a - b).abs() <= 1e-10_f64.max(1e-12 * a.abs().max(b.abs())),
                "{context}: {a} != {b}"
            );
        }
        (Value::Array(a), Value::Array(b)) => {
            assert_eq!(a.len(), b.len(), "{context}");
            for (i, (a, b)) in a.iter().zip(b).enumerate() {
                compare(a, b, &format!("{context}[{i}]"));
            }
        }
        (Value::Object(a), Value::Object(b)) => {
            assert_eq!(
                a.keys().collect::<Vec<_>>(),
                b.keys().collect::<Vec<_>>(),
                "{context}"
            );
            for (key, a) in a {
                compare(a, &b[key], &format!("{context}.{key}"));
            }
        }
        _ => assert_eq!(actual, expected, "{context}"),
    }
}

#[test]
fn every_rule_and_setting_matches_independent_toolkit_interpretation() {
    use sha2::{Digest, Sha256};
    let report = fixture();
    assert_eq!(
        format!(
            "{:x}",
            Sha256::digest(include_bytes!("../../data/rosemary.offxml"))
        ),
        report["source_sha256"]
    );
    assert_eq!(report["toolkit_version"], "0.19.0");
    assert_eq!(report["records"].as_array().unwrap().len(), 6);
    for record in report["records"].as_array().unwrap() {
        let ff = ForceField::from_offxml(record["xml"].as_str().unwrap()).unwrap();
        compare(
            &snapshot(&ff),
            &json!({"parameters":record["parameters"],"settings":record["settings"]}),
            record["name"].as_str().unwrap(),
        );
    }
}

#[test]
fn optional_sections_are_empty_without_changing_other_rules() {
    let report = fixture();
    let record = report["records"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["name"] == "optional-sections-absent")
        .unwrap();
    let ff = ForceField::from_offxml(record["xml"].as_str().unwrap()).unwrap();
    assert!(ff.constraints.is_empty() && ff.impropers.is_empty() && ff.library.is_empty());
    let reference = crate::reference_records();
    for record in reference {
        let input = kekule::smiles::to_molecules(record["mapped_smiles"].as_str().unwrap())
            .unwrap()
            .remove(0);
        let actual = ff.label_molecule(&input).unwrap();
        let original = ForceField::rosemary()
            .unwrap()
            .label_molecule(&input)
            .unwrap();
        assert!(actual.constraints.is_empty() && actual.improper_torsions.is_empty());
        assert_eq!(actual.bonds, original.bonds);
        assert_eq!(actual.angles, original.angles);
        assert_eq!(actual.proper_torsions, original.proper_torsions);
        assert_eq!(actual.vdw, original.vdw);
    }
}
