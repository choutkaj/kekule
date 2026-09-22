//! Optional validation consumer, not a production dynamics backend.
//! nm, kJ/mol, e; vacuum NoCutoff; retain constrained valence forces.
use kekule::topology::InstanceAtomId;
use kekule_openff::ParameterizedTopology;
use potentials::base::{Potential2, Potential3, Potential4};
use serde_json::{json, Value};
use std::collections::BTreeMap;

// OpenMM's electrostatic conversion, kJ mol^-1 nm e^-2.
const COULOMB: f64 = 138.935_457_644_381_98;
type Xyz = [f64; 3];
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
    let axis = axis.map(|v| v / norm);
    let v = sub(a, b);
    let w = sub(d, c);
    let v = sub(v, axis.map(|x| x * dot(v, axis)));
    let w = sub(w, axis.map(|x| x * dot(w, axis)));
    let scale = (dot(v, v) * dot(w, w)).sqrt();
    if scale < 1e-18 {
        return Err("degenerate torsion plane".into());
    }
    Ok((
        (dot(v, w) / scale).clamp(-1.0, 1.0),
        (dot(cross(axis, v), w) / scale).clamp(-1.0, 1.0),
    ))
}

fn evaluate(system: &ParameterizedTopology, xyz: &[Xyz], charges: &[f64]) -> Result<Value, String> {
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
        bonds += potentials::bond::Harm::<f64>::new(
            *term.parameter.k.value() / 2.0,
            *term.parameter.length.value(),
        )
        .energy(distance(a, b));
    }
    let mut angles = 0.0;
    for term in system.angles() {
        let [a, b, c] = term.atoms.map(|a| index[&a]);
        let u = sub(xyz[a], xyz[b]);
        let v = sub(xyz[c], xyz[b]);
        let r1 = dot(u, u);
        let r2 = dot(v, v);
        if r1 * r2 < 1e-24 {
            return Err("degenerate angle".into());
        }
        let cosine = (dot(u, v) / (r1 * r2).sqrt()).clamp(-1.0, 1.0);
        angles += potentials::angle::Harm::<f64>::new(
            *term.parameter.k.value() / 2.0,
            *term.parameter.angle.value(),
        )
        .energy(r1, r2, cosine);
    }
    let mut torsions = [0.0; 2];
    for (sum, terms) in torsions
        .iter_mut()
        .zip([system.proper_torsions(), system.improper_torsions()])
    {
        for term in terms {
            let (cosine, sine) = dihedral(term.atoms.map(|a| xyz[index[&a]]))?;
            for p in &term.parameter.terms {
                *sum += potentials::torsion::Cos::<f64>::new(
                    *p.k.value() / p.idivf,
                    p.periodicity as i32,
                    *p.phase.value(),
                )
                .energy(cosine, sine);
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
                vdw += vs * potentials::pair::Lj::<f64>::new(epsilon, sigma).energy(r2);
            }
            if qs != 0.0 {
                electrostatics +=
                    potentials::pair::Coul::<f64>::new(COULOMB * charges[i] * charges[j] * qs)
                        .energy(r2);
            }
        }
    }
    let components = [bonds, angles, torsions[0], torsions[1], vdw, electrostatics];
    if components.iter().any(|v| !v.is_finite()) {
        return Err("nonfinite energy".into());
    }
    Ok(
        json!({"Bonds":bonds,"Angles":angles,"ProperTorsions":torsions[0],"ImproperTorsions":torsions[1],"vdW":vdw,"Electrostatics":electrostatics,"Total":components.iter().sum::<f64>()}),
    )
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
        result.push(json!({"native_charges":evaluate(system,&xyz,system.charges().value())?,"reference_charges":evaluate(system,&xyz,&reference)?}));
    }
    Ok(json!(result))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn signed_dihedral_and_phase_are_preserved() {
        let (c, s) = dihedral([[1., 0., 0.], [0., 0., 0.], [0., 0., 1.], [0., 1., 1.]]).unwrap();
        assert!(c.abs() < 1e-14 && (s - 1.).abs() < 1e-14);
        let p = potentials::torsion::Cos::<f64>::new(2., 1, std::f64::consts::FRAC_PI_2);
        assert!((p.energy(c, s) - 4.).abs() < 1e-14);
        assert!(dihedral([[0.; 3]; 4]).is_err());
    }
}
