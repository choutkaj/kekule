//! Native observer for the optional externally supplied OpenFF reference corpus.
use kekule::{
    core::{AtomId, Molecule},
    smiles,
    units::Quantity,
};
use kekule_openff::{ForceField, Interaction, NaglModel, ParameterIdentity};
use serde_json::{json, Value};
use std::io::{self, BufRead};
#[path = "openff_parameterize/energy.rs"]
mod energy;

fn labels<const N: usize>(
    map: &std::collections::BTreeMap<[AtomId; N], ParameterIdentity>,
    m: &Molecule,
) -> Vec<Value> {
    map.iter().map(|(a,p)|json!({"maps":a.iter().map(|a|m.atom(*a).unwrap().atom_map.unwrap()).collect::<Vec<_>>(),"id":p.id,"smirks":p.smirks})).collect()
}
fn terms<const N: usize, P>(
    terms: &[Interaction<N, P>],
    m: &Molecule,
    value: impl Fn(&P) -> Value,
) -> Vec<Value> {
    terms.iter().map(|t|json!({"maps":t.atoms.iter().map(|a|m.atom(a.atom()).unwrap().atom_map.unwrap()).collect::<Vec<_>>(),"parameter":value(&t.parameter)})).collect()
}
fn observe(
    ff: &ForceField,
    model: Option<&NaglModel>,
    source: &str,
    mode: &str,
    request: &Value,
) -> Result<Value, Box<dyn std::error::Error>> {
    let mut molecules = smiles::to_molecules(source)?;
    if molecules.len() != 1 {
        return Err("expected one molecule".into());
    }
    let molecule = molecules.remove(0);
    let maps = molecule
        .atoms()
        .map(|(_, a)| a.atom_map.ok_or("missing atom map"))
        .collect::<Result<Vec<_>, _>>()?;
    let labels_ = ff.label_molecule(&molecule)?;
    let mut result = json!({"status":"ok","maps":maps,"labels":{
        "Bonds":labels(&labels_.bonds,&molecule),"Angles":labels(&labels_.angles,&molecule),
        "ProperTorsions":labels(&labels_.proper_torsions,&molecule),"ImproperTorsions":labels(&labels_.improper_torsions,&molecule),
        "Constraints":labels(&labels_.constraints,&molecule),"vdW":labels(&labels_.vdw,&molecule)}});
    if let Some(model) = model {
        match kekule_openff::diagnostics::lookup_key(model, &molecule) {
            Ok(key) => result["lookup_key"] = json!(key),
            Err(e) => {
                result["lookup_key"] = Value::Null;
                result["identity_error"] = json!(e.to_string());
            }
        }
        result["features"] = match kekule_openff::diagnostics::atom_features(model, &molecule) {
            Ok(f) => json!({"status":"ok","values":f}),
            Err(e) => json!({"status":"error","message":e.to_string()}),
        };
        if mode != "validate" {
            result["charges"] = match model.assign_charges(&molecule) {
                Ok(q) => {
                    json!({"status":"ok","values":q.charges.value(),"source":format!("{:?}",q.source)})
                }
                Err(e) => json!({"status":"error","message":e.to_string()}),
            };
            result["inference"] = match kekule_openff::diagnostics::infer_charges(model, &molecule)
            {
                Ok(q) => json!({"status":"ok","values":q.charges.value()}),
                Err(e) => json!({"status":"error","message":e.to_string()}),
            };
        }
        let system = ff.parameterize_molecule(molecule.clone(), model)?;
        if !request["coordinates_nm"].is_null() {
            result["energies"] = energy::observe(&system, &maps, request)?;
        }
        result["charge_sources"] = json!(system
            .charge_sources()
            .iter()
            .map(|s| format!("{s:?}"))
            .collect::<Vec<_>>());
        let scalar = |q: &Quantity<f64>| json!(q.value());
        let torsion = |p: &kekule_openff::TorsionParameter| json!({"id":p.source.id,"terms":p.terms.iter().map(|t|json!({"k":scalar(&t.k),"phase":scalar(&t.phase),"periodicity":t.periodicity,"idivf":t.idivf})).collect::<Vec<_>>()});
        result["system"] = json!({"charges":system.charges().value(),
            "exceptions":system.pair_exceptions().iter().map(|p|json!({"maps":p.atoms.map(|a|molecule.atom(a.atom()).unwrap().atom_map.unwrap()),"vdw_scale":p.vdw_scale,"electrostatics_scale":p.electrostatics_scale})).collect::<Vec<_>>(),
            "bonds":terms(system.bonds(),&molecule,|p|json!({"id":p.source.id,"length":scalar(&p.length),"k":scalar(&p.k)})),
            "angles":terms(system.angles(),&molecule,|p|json!({"id":p.source.id,"angle":scalar(&p.angle),"k":scalar(&p.k)})),
            "constraints":terms(system.constraints(),&molecule,|p|json!({"id":p.source.id,"distance":scalar(&p.distance)})),
            "propers":terms(system.proper_torsions(),&molecule,torsion),"impropers":terms(system.improper_torsions(),&molecule,torsion),
            "vdw":system.vdw().iter().map(|p|json!({"id":p.source.id,"sigma":scalar(&p.sigma),"epsilon":scalar(&p.epsilon)})).collect::<Vec<_>>()});
    }
    Ok(result)
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let ff = std::env::args()
        .nth(2)
        .map(ForceField::from_file)
        .unwrap_or_else(ForceField::rosemary)?;
    let model = std::env::args().nth(1).map(NaglModel::load).transpose()?;
    for line in io::stdin().lock().lines() {
        let line = line?;
        let result = serde_json::from_str::<Value>(&line)
            .map_err(|e| e.to_string())
            .and_then(|v| {
                observe(
                    &ff,
                    model.as_ref(),
                    v["smiles"].as_str().ok_or("missing smiles")?,
                    v["mode"].as_str().unwrap_or("full"),
                    &v,
                )
                .map_err(|e| e.to_string())
            });
        println!(
            "{}",
            result.unwrap_or_else(|e| json!({"status":"error","message":e}))
        );
    }
    Ok(())
}
