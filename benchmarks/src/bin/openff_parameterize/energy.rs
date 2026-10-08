//! Optional validation consumer of the public `kekule-potentials` OpenFF
//! evaluator: nm, kJ/mol, e; vacuum NoCutoff; constrained bond terms retained.
//! Cartesian gradients are dE/dx in kJ/mol/nm, not forces. Near-linear
//! derivatives are not capped to mimic a reference backend.
use kekule::{
    geometry::Point3,
    structure::{Model, Positions},
    units::{Quantity, ELEMENTARY_CHARGE, NANOMETER},
};
use kekule_openff::ParameterizedTopology;
use kekule_potentials::{openff::OpenFfPotential, ComponentKind, Potential};
use serde_json::{json, Value};
use std::sync::Arc;

type Xyz = [f64; 3];
const COMPONENTS: [(ComponentKind, &str); 6] = [
    (ComponentKind::Bonds, "Bonds"),
    (ComponentKind::Angles, "Angles"),
    (ComponentKind::ProperTorsions, "ProperTorsions"),
    (ComponentKind::ImproperTorsions, "ImproperTorsions"),
    (ComponentKind::VanDerWaals, "vdW"),
    (ComponentKind::Electrostatics, "Electrostatics"),
];

fn evaluate(
    potential: &OpenFfPotential,
    system: &ParameterizedTopology,
    xyz: &[Xyz],
    derivatives: bool,
) -> Result<Value, String> {
    let points: Vec<Point3> = xyz.iter().map(|&[x, y, z]| Point3::new(x, y, z)).collect();
    let positions =
        Positions::from_vec(Quantity::new(points, NANOMETER)).map_err(|e| e.to_string())?;
    let model = Model::new(Arc::clone(system.topology()), positions).map_err(|e| e.to_string())?;
    let energy = if derivatives {
        potential.evaluate(model.view()).map(|e| e.energy().clone())
    } else {
        potential.energy(model.view())
    }
    .map_err(|e| e.to_string())?;
    let mut result = serde_json::Map::new();
    for (kind, name) in COMPONENTS {
        let value = energy
            .component(kind)
            .ok_or_else(|| format!("missing {kind} component"))?;
        result.insert(name.into(), json!(value.value()));
    }
    result.insert("Total".into(), json!(energy.total().value()));
    let mut result = Value::Object(result);
    if derivatives {
        let vectors = |g: &[kekule::geometry::Vector3]| {
            json!(g.iter().map(|v| [v.x, v.y, v.z]).collect::<Vec<_>>())
        };
        let mut values = serde_json::Map::new();
        for (component, (kind, name)) in potential
            .evaluate_components(model.view())
            .map_err(|e| e.to_string())?
            .into_iter()
            .zip(COMPONENTS)
        {
            if component.kind != kind {
                return Err(format!("unexpected component order at {name}"));
            }
            values.insert(name.into(), vectors(component.gradient.value()));
        }
        let total = potential
            .evaluate(model.view())
            .map_err(|e| e.to_string())?;
        values.insert("Total".into(), vectors(total.gradient().value()));
        result["gradients"] = Value::Object(values);
    }
    Ok(result)
}

pub(super) fn observe(
    system: &ParameterizedTopology,
    maps: &[u32],
    request: &Value,
) -> Result<Value, String> {
    let frames: Vec<Vec<Xyz>> =
        serde_json::from_value(request["coordinates_nm"].clone()).map_err(|e| e.to_string())?;
    let reference: Vec<f64> =
        serde_json::from_value(request["reference_charges"].clone()).map_err(|e| e.to_string())?;
    let n = maps.len();
    let mut sorted = maps.to_vec();
    sorted.sort();
    if sorted != (1..=n as u32).collect::<Vec<_>>()
        || reference.len() != n
        || reference.iter().any(|v| !v.is_finite())
    {
        return Err(
            "energy observation requires consecutive unique maps and finite reference charges"
                .into(),
        );
    }
    let reference: Vec<_> = maps.iter().map(|&m| reference[m as usize - 1]).collect();
    let native = OpenFfPotential::new(system).map_err(|e| e.to_string())?;
    let referenced = native
        .clone()
        .with_charges(Quantity::new(reference, ELEMENTARY_CHARGE))
        .map_err(|e| e.to_string())?;
    let mut result = vec![];
    for frame in frames {
        if frame.len() != n || frame.iter().flatten().any(|x| !x.is_finite()) {
            return Err("invalid coordinates".into());
        }
        let xyz: Vec<_> = maps.iter().map(|&m| frame[m as usize - 1]).collect();
        let derivatives = request["gradients"].as_bool().unwrap_or(false);
        let mut observed = json!({"native_charges":evaluate(&native,system,&xyz,derivatives)?,"reference_charges":evaluate(&referenced,system,&xyz,derivatives)?});
        if let Some(directions) = request["directions"].as_array() {
            let mut differences = vec![];
            for direction in directions {
                let direction: Vec<Xyz> =
                    serde_json::from_value(direction.clone()).map_err(|e| e.to_string())?;
                if direction.len() != n || direction.iter().flatten().any(|x| !x.is_finite()) {
                    return Err("invalid finite difference direction".into());
                }
                let direction: Vec<_> = maps.iter().map(|&m| direction[m as usize - 1]).collect();
                let mut steps = vec![];
                let step_sizes: Vec<f64> = if request["finite_difference_steps_nm"].is_null() {
                    vec![1e-5, 5e-6]
                } else {
                    serde_json::from_value(request["finite_difference_steps_nm"].clone())
                        .map_err(|e| e.to_string())?
                };
                if step_sizes.is_empty() || step_sizes.iter().any(|h| !h.is_finite() || *h <= 0.0) {
                    return Err("invalid finite difference steps".into());
                }
                for h in step_sizes {
                    let mut energies = vec![];
                    for sign in [-1.0, 1.0] {
                        let shifted: Vec<Xyz> = xyz
                            .iter()
                            .zip(&direction)
                            .map(|(x, d)| std::array::from_fn(|j| x[j] + sign * h * d[j]))
                            .collect();
                        energies.push(evaluate(&referenced, system, &shifted, false)?);
                    }
                    steps.push(json!({"step_nm":h,"minus":energies[0],"plus":energies[1]}));
                }
                differences.push(json!(steps));
            }
            observed["finite_differences"] = json!(differences);
        }
        result.push(observed);
    }
    Ok(json!(result))
}
