//! Optional validation consumer, not a production dynamics backend.
//! nm, kJ/mol, e; vacuum NoCutoff; retain constrained valence forces.
//! Optional Cartesian gradients are dE/dx in kJ/mol/nm, not forces.
//! Near-linear derivatives are not capped to mimic a reference backend.
use kekule::topology::InstanceAtomId;
use kekule_openff::ParameterizedTopology;
use potentials::base::{Potential2, Potential4};
use serde_json::{json, Value};
use std::collections::BTreeMap;

// OpenMM's electrostatic conversion, kJ mol^-1 nm e^-2.
const COULOMB: f64 = 138.935_457_644_381_98;
type Xyz = [f64; 3];
const COMPONENTS: [&str; 6] = [
    "Bonds",
    "Angles",
    "ProperTorsions",
    "ImproperTorsions",
    "vdW",
    "Electrostatics",
];

fn add(gradient: &mut [Xyz], atom: usize, vector: Xyz, scale: f64) {
    for (g, v) in gradient[atom].iter_mut().zip(vector) {
        *g += scale * v;
    }
}
fn sub(a: Xyz, b: Xyz) -> Xyz {
    std::array::from_fn(|i| a[i] - b[i])
}
fn dot(a: Xyz, b: Xyz) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: Xyz, b: Xyz) -> Xyz {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dihedral([a, b, c, d]: [Xyz; 4]) -> Result<(f64, f64), String> {
    let axis = sub(c, b);
    let norm = dot(axis, axis).sqrt();
    if norm < 1e-12 {
        return Err("degenerate torsion axis".into());
    }
    let v = sub(a, b);
    let w = sub(d, c);
    let v = cross(axis, v);
    let w = cross(axis, w);
    let scale = (dot(v, v) * dot(w, w)).sqrt();
    if scale < 1e-18 * norm * norm {
        return Err("degenerate torsion plane".into());
    }
    Ok((
        (dot(v, w) / scale).clamp(-1.0, 1.0),
        (dot(axis, cross(v, w)) / (norm * scale)).clamp(-1.0, 1.0),
    ))
}

fn dihedral_gradient([a, b, c, d]: [Xyz; 4]) -> [Xyz; 4] {
    let left = sub(a, b);
    let axis = sub(c, b);
    let right = sub(d, c);
    let axis_sq = dot(axis, axis);
    let n0 = cross(left, axis);
    let n3 = cross(axis, right);
    let g0 = n0.map(|v| v * axis_sq.sqrt() / dot(n0, n0));
    let g3 = n3.map(|v| v * axis_sq.sqrt() / dot(n3, n3));
    let alpha = dot(left, axis) / axis_sq;
    let beta = dot(right, axis) / axis_sq;
    [
        g0,
        std::array::from_fn(|i| (alpha - 1.0) * g0[i] + beta * g3[i]),
        std::array::from_fn(|i| -alpha * g0[i] - (beta + 1.0) * g3[i]),
        g3,
    ]
}

fn angle([a, b, c]: [Xyz; 3], k: f64, theta0: f64) -> Result<(f64, [Xyz; 3]), String> {
    let u = sub(a, b);
    let v = sub(c, b);
    let r1 = dot(u, u);
    let r2 = dot(v, v);
    if r1 * r2 < 1e-24 {
        return Err("degenerate angle".into());
    }
    let normal = cross(u, v);
    let norm = dot(normal, normal).sqrt();
    // atan2 avoids acos roundoff near linearity. The harmonic energy is
    // unchanged mathematically. Neither its gradient nor its normal is capped.
    let delta = norm.atan2(dot(u, v)) - theta0;
    let energy = 0.5 * k * delta * delta;
    if norm == 0.0 {
        if delta.abs() > 1e-12 && k != 0.0 {
            return Err("undefined gradient at a nonstationary linear angle".into());
        }
        return Ok((energy, [[0.; 3]; 3]));
    }
    let normal = normal.map(|x| x / norm);
    let da = cross(u, normal).map(|x| k * delta * x / r1);
    let dc = cross(normal, v).map(|x| k * delta * x / r2);
    Ok((energy, [da, std::array::from_fn(|i| -da[i] - dc[i]), dc]))
}

fn evaluate(
    system: &ParameterizedTopology,
    xyz: &[Xyz],
    charges: &[f64],
    derivatives: bool,
) -> Result<Value, String> {
    let mut gradients: [Vec<Xyz>; 6] =
        std::array::from_fn(|_| vec![[0.0; 3]; if derivatives { xyz.len() } else { 0 }]);
    let index: BTreeMap<InstanceAtomId, usize> = system
        .topology()
        .atom_ids()
        .iter()
        .copied()
        .enumerate()
        .map(|(i, a)| (a, i))
        .collect();
    let distance = |a: usize, b: usize| {
        let r = sub(xyz[a], xyz[b]);
        dot(r, r)
    };
    let mut bonds = 0.0;
    for term in system.bonds() {
        let [a, b] = term.atoms.map(|a| index[&a]);
        let p = potentials::bond::Harm::<f64>::new(
            *term.parameter.k.value() / 2.0,
            *term.parameter.length.value(),
        );
        let r2 = distance(a, b);
        bonds += p.energy(r2);
        if derivatives {
            if r2 < 1e-24 {
                return Err("coincident bonded atoms".into());
            }
            let force = p.force_factor(r2);
            add(&mut gradients[0], a, sub(xyz[a], xyz[b]), -force);
            add(&mut gradients[0], b, sub(xyz[a], xyz[b]), force);
        }
    }
    let mut angles = 0.0;
    for term in system.angles() {
        let [a, b, c] = term.atoms.map(|a| index[&a]);
        let (energy, gradient) = angle(
            [xyz[a], xyz[b], xyz[c]],
            *term.parameter.k.value(),
            *term.parameter.angle.value(),
        )?;
        angles += energy;
        if derivatives {
            for (atom, vector) in [a, b, c].into_iter().zip(gradient) {
                add(&mut gradients[1], atom, vector, 1.0);
            }
        }
    }
    let mut torsions = [0.0; 2];
    for (component, (sum, terms)) in torsions
        .iter_mut()
        .zip([system.proper_torsions(), system.improper_torsions()])
        .enumerate()
    {
        for term in terms {
            let atoms = term.atoms.map(|a| index[&a]);
            let points = atoms.map(|a| xyz[a]);
            let (cosine, sine) = dihedral(points)?;
            for p in &term.parameter.terms {
                let p = potentials::torsion::Cos::<f64>::new(
                    *p.k.value() / p.idivf,
                    p.periodicity as i32,
                    *p.phase.value(),
                );
                *sum += p.energy(cosine, sine);
                if derivatives {
                    let derivative = p.derivative(cosine, sine);
                    for (atom, vector) in atoms.into_iter().zip(dihedral_gradient(points)) {
                        add(&mut gradients[component + 2], atom, vector, derivative);
                    }
                }
            }
        }
    }
    let exceptions: BTreeMap<_, _> = system
        .pair_exceptions()
        .iter()
        .map(|p| {
            let mut atoms = p.atoms.map(|a| index[&a]);
            atoms.sort();
            (atoms, (p.vdw_scale, p.electrostatics_scale))
        })
        .collect();
    let (mut vdw, mut electrostatics) = (0.0, 0.0);
    for i in 0..xyz.len() {
        for j in i + 1..xyz.len() {
            let (vs, qs) = exceptions.get(&[i, j]).copied().unwrap_or((1.0, 1.0));
            if vs == 0.0 && qs == 0.0 {
                continue;
            }
            let r2 = distance(i, j);
            if r2 < 1e-24 {
                return Err("coincident nonbonded atoms".into());
            }
            let a = &system.vdw()[i];
            let b = &system.vdw()[j];
            if vs != 0.0 {
                let sigma = (a.sigma.value() + b.sigma.value()) / 2.0;
                let epsilon = (a.epsilon.value() * b.epsilon.value()).sqrt();
                let p = potentials::pair::Lj::<f64>::new(epsilon, sigma);
                vdw += vs * p.energy(r2);
                if derivatives {
                    let force = vs * p.force_factor(r2);
                    add(&mut gradients[4], i, sub(xyz[i], xyz[j]), -force);
                    add(&mut gradients[4], j, sub(xyz[i], xyz[j]), force);
                }
            }
            if qs != 0.0 {
                let p = potentials::pair::Coul::<f64>::new(COULOMB * charges[i] * charges[j] * qs);
                electrostatics += p.energy(r2);
                if derivatives {
                    let force = p.force_factor(r2);
                    add(&mut gradients[5], i, sub(xyz[i], xyz[j]), -force);
                    add(&mut gradients[5], j, sub(xyz[i], xyz[j]), force);
                }
            }
        }
    }
    let components = [bonds, angles, torsions[0], torsions[1], vdw, electrostatics];
    if components.iter().any(|v| !v.is_finite()) {
        return Err("nonfinite energy".into());
    }
    let mut result = json!({"Bonds":bonds,"Angles":angles,"ProperTorsions":torsions[0],"ImproperTorsions":torsions[1],"vdW":vdw,"Electrostatics":electrostatics,"Total":components.iter().sum::<f64>()});
    if derivatives {
        if gradients.iter().flatten().flatten().any(|v| !v.is_finite()) {
            return Err("nonfinite gradient".into());
        }
        let total: Vec<Xyz> = (0..xyz.len())
            .map(|i| std::array::from_fn(|j| gradients.iter().map(|g| g[i][j]).sum::<f64>()))
            .collect();
        let mut values = serde_json::Map::new();
        for (name, gradient) in COMPONENTS.into_iter().zip(gradients) {
            values.insert(name.into(), json!(gradient));
        }
        values.insert("Total".into(), json!(total));
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
    let mut result = vec![];
    for frame in frames {
        if frame.len() != n || frame.iter().flatten().any(|x| !x.is_finite()) {
            return Err("invalid coordinates".into());
        }
        let xyz: Vec<_> = maps.iter().map(|&m| frame[m as usize - 1]).collect();
        let derivatives = request["gradients"].as_bool().unwrap_or(false);
        let mut observed = json!({"native_charges":evaluate(system,&xyz,system.charges().value(),derivatives)?,"reference_charges":evaluate(system,&xyz,&reference,derivatives)?});
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
                        energies.push(evaluate(system, &shifted, &reference, false)?);
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn externally_sourced_near_linear_angle_differentiates_its_energy() {
        // PubChem CID 443915, unchanged geometry from the full validation panel.
        // This catches both the sin(theta) floor and acos cancellation.
        let bytes = std::fs::read(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/openff/data/reference.json.gz"
        ))
        .unwrap();
        let reference: Value =
            serde_json::from_reader(flate2::read::GzDecoder::new(&bytes[..])).unwrap();
        let case = reference["records"]
            .as_array()
            .unwrap()
            .iter()
            .find(|r| r["id"] == "pubchem-443915")
            .unwrap();
        let parameter = case["parameters"]["Angles"]
            .as_array()
            .unwrap()
            .iter()
            .find(|p| p["atoms"] == json!([16, 17, 18]))
            .unwrap();
        let points: [Xyz; 3] = std::array::from_fn(|i| {
            serde_json::from_value(case["coordinates_nm"][0][15 + i].clone()).unwrap()
        });
        let k = parameter["values"]["k"].as_f64().unwrap();
        let theta0 = parameter["values"]["angle"].as_f64().unwrap();
        let (_, gradient) = angle(points, k, theta0).unwrap();
        for atom in 0..3 {
            for axis in 0..3 {
                let estimates: Vec<_> = [1e-7, 5e-8]
                    .into_iter()
                    .map(|h| {
                        let mut plus = points;
                        let mut minus = points;
                        plus[atom][axis] += h;
                        minus[atom][axis] -= h;
                        (angle(plus, k, theta0).unwrap().0 - angle(minus, k, theta0).unwrap().0)
                            / (2.0 * h)
                    })
                    .collect();
                let numerical = (4.0 * estimates[1] - estimates[0]) / 3.0;
                assert!(
                    (gradient[atom][axis] - numerical).abs() < 2e-4 + 1e-7 * numerical.abs(),
                    "{atom}/{axis}: analytic={} numerical={numerical}",
                    gradient[atom][axis]
                );
            }
        }
        assert!(angle([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]], k, theta0).is_err());
    }

    #[test]
    fn signed_dihedral_and_phase_are_preserved() {
        let (c, s) = dihedral([[1., 0., 0.], [0., 0., 0.], [0., 0., 1.], [0., 1., 1.]]).unwrap();
        assert!(c.abs() < 1e-14 && (s - 1.).abs() < 1e-14);
        let p = potentials::torsion::Cos::<f64>::new(2., 1, std::f64::consts::FRAC_PI_2);
        assert!((p.energy(c, s) - 4.).abs() < 1e-14);
        assert!(dihedral([[0.; 3]; 4]).is_err());
    }
}
