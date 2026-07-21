use dreid_kernel::potentials::bonded::{
    CosineHarmonic, CosineLinear, Harmonic, PlanarInversion, Torsion, UmbrellaInversion,
};
use dreid_kernel::potentials::nonbonded::{Coulomb, HydrogenBond, LennardJones};
use dreid_kernel::{AngleKernel, HybridKernel, PairKernel, TorsionKernel};
use molecular::core::Point3;
use molecular::modeling::potential::{
    Potential, PotentialError, PotentialEvaluation, PotentialGeometryError, Vector3,
};
use molecular::modeling::{InstanceAtomId, Model};
use molecular::units::{Quantity, MODEL_ENERGY_UNIT, MODEL_GRADIENT_UNIT};

use crate::geometry::{
    angle_cosine, displacement, hydrogen_bond_cosine, inversion_cosine, torsion, GeometryError,
};
use crate::prepare::{AngleTerm, DreidingPotential, InversionTerm};

impl Potential for DreidingPotential {
    fn evaluate(&mut self, model: &Model) -> Result<PotentialEvaluation, PotentialError> {
        if &self.definition != model.definition_key() {
            return Err(PotentialError::IncompatibleModel);
        }

        let positions = model.positions().into_value();
        let mut energy = 0.0;
        let mut gradient = vec![Vector3::zero(); positions.len()];

        for term in &self.bonds {
            let (vector, distance_sq) = displacement(positions[term.a], positions[term.b])
                .map_err(|error| {
                    geometry_error(&self.atom_ids, "DREIDING bond", &[term.a, term.b], error)
                })?;
            let result = Harmonic::compute(distance_sq, (term.k_half, term.r0));
            energy += result.energy;
            add_pair_gradient(&mut gradient, term.a, term.b, vector, result.diff);
        }

        for term in &self.angles {
            let (atoms, result, derivative) = match *term {
                AngleTerm::Harmonic {
                    atoms,
                    c_half,
                    cos0,
                } => {
                    let (cosine, derivative) = angle_geometry(&self.atom_ids, positions, atoms)?;
                    (
                        atoms,
                        CosineHarmonic::compute(cosine, (c_half, cos0)),
                        derivative,
                    )
                }
                AngleTerm::Linear { atoms, c } => {
                    let (cosine, derivative) = angle_geometry(&self.atom_ids, positions, atoms)?;
                    (atoms, CosineLinear::compute(cosine, c), derivative)
                }
            };
            energy += result.energy;
            add_derivative::<9, 3>(&mut gradient, atoms, derivative, result.diff);
        }

        for term in &self.torsions {
            let points = points4(positions, term.atoms);
            let (cosine, sine, derivative) = torsion(points).map_err(|error| {
                geometry_error(&self.atom_ids, "DREIDING torsion", &term.atoms, error)
            })?;
            let result = Torsion::compute(
                cosine,
                sine,
                (term.v_half, term.n, term.cos_n_phi0, term.sin_n_phi0),
            );
            energy += result.energy;
            add_derivative::<12, 4>(&mut gradient, term.atoms, derivative, result.diff);
        }

        for term in &self.inversions {
            let (atoms, result, derivative) = match *term {
                InversionTerm::Planar { atoms, c_half } => {
                    let (cosine, derivative) = inversion_cosine(points4(positions, atoms))
                        .map_err(|error| {
                            geometry_error(&self.atom_ids, "DREIDING inversion", &atoms, error)
                        })?;
                    (atoms, PlanarInversion::compute(cosine, c_half), derivative)
                }
                InversionTerm::Umbrella {
                    atoms,
                    c_half,
                    cos_psi0,
                } => {
                    let (cosine, derivative) = inversion_cosine(points4(positions, atoms))
                        .map_err(|error| {
                            geometry_error(&self.atom_ids, "DREIDING inversion", &atoms, error)
                        })?;
                    (
                        atoms,
                        UmbrellaInversion::compute(cosine, (c_half, cos_psi0)),
                        derivative,
                    )
                }
            };
            energy += result.energy;
            add_derivative::<12, 4>(&mut gradient, atoms, derivative, result.diff);
        }

        for term in &self.nonbonded {
            let (vector, distance_sq) = displacement(positions[term.first], positions[term.second])
                .map_err(|error| {
                    geometry_error(
                        &self.atom_ids,
                        "DREIDING nonbonded pair",
                        &[term.first, term.second],
                        error,
                    )
                })?;
            let vdw = LennardJones::compute(distance_sq, (term.d0, term.r0_sq));
            let electrostatic = Coulomb::compute(distance_sq, term.coulomb);
            energy += vdw.energy + electrostatic.energy;
            add_pair_gradient(
                &mut gradient,
                term.first,
                term.second,
                vector,
                vdw.diff + electrostatic.diff,
            );
        }

        for term in &self.hydrogen_bonds {
            let donor = positions[term.donor];
            let acceptor = positions[term.acceptor];
            let (vector, distance_sq) = displacement(donor, acceptor).map_err(|error| {
                geometry_error(
                    &self.atom_ids,
                    "DREIDING hydrogen bond",
                    &[term.donor, term.hydrogen, term.acceptor],
                    error,
                )
            })?;
            let (cosine, derivative) =
                hydrogen_bond_cosine(donor, positions[term.hydrogen], acceptor).map_err(
                    |error| {
                        geometry_error(
                            &self.atom_ids,
                            "DREIDING hydrogen bond",
                            &[term.donor, term.hydrogen, term.acceptor],
                            error,
                        )
                    },
                )?;
            let result = HydrogenBond::<4>::compute(distance_sq, cosine, (term.d_hb, term.r_hb_sq));
            energy += result.energy;
            add_pair_gradient(
                &mut gradient,
                term.donor,
                term.acceptor,
                vector,
                result.force_factor_rad,
            );
            add_derivative::<9, 3>(
                &mut gradient,
                [term.donor, term.hydrogen, term.acceptor],
                derivative,
                result.force_factor_ang,
            );
        }

        PotentialEvaluation::new(
            model,
            Quantity::new(energy, MODEL_ENERGY_UNIT),
            Quantity::new(gradient, MODEL_GRADIENT_UNIT),
        )
    }
}

fn angle_geometry(
    atom_ids: &[InstanceAtomId],
    positions: &[Point3],
    atoms: [usize; 3],
) -> Result<(f64, [f64; 9]), PotentialError> {
    angle_cosine(
        positions[atoms[0]],
        positions[atoms[1]],
        positions[atoms[2]],
    )
    .map_err(|error| geometry_error(atom_ids, "DREIDING angle", &atoms, error))
}

fn points4(positions: &[Point3], atoms: [usize; 4]) -> [Point3; 4] {
    [
        positions[atoms[0]],
        positions[atoms[1]],
        positions[atoms[2]],
        positions[atoms[3]],
    ]
}

fn add_pair_gradient(
    gradient: &mut [Vector3],
    first: usize,
    second: usize,
    displacement: [f64; 3],
    force_factor: f64,
) {
    add_vector(&mut gradient[first], displacement, -force_factor);
    add_vector(&mut gradient[second], displacement, force_factor);
}

fn add_derivative<const N: usize, const A: usize>(
    gradient: &mut [Vector3],
    atoms: [usize; A],
    derivative: [f64; N],
    factor: f64,
) {
    debug_assert_eq!(N, A * 3);
    for (slot, atom) in atoms.into_iter().enumerate() {
        gradient[atom].x += factor * derivative[slot * 3];
        gradient[atom].y += factor * derivative[slot * 3 + 1];
        gradient[atom].z += factor * derivative[slot * 3 + 2];
    }
}

fn add_vector(target: &mut Vector3, vector: [f64; 3], scale: f64) {
    target.x += vector[0] * scale;
    target.y += vector[1] * scale;
    target.z += vector[2] * scale;
}

fn geometry_error(
    atom_ids: &[InstanceAtomId],
    interaction: &'static str,
    atoms: &[usize],
    error: GeometryError,
) -> PotentialError {
    let kind = match error {
        GeometryError::Coincident => PotentialGeometryError::CoincidentAtoms,
        GeometryError::DegenerateAngle => PotentialGeometryError::DegenerateAngle,
        GeometryError::DegenerateDihedral => PotentialGeometryError::DegenerateDihedral,
        GeometryError::DegenerateInversion => PotentialGeometryError::DegenerateInversion,
    };
    let atoms = atoms.iter().map(|&atom| atom_ids[atom]).collect::<Vec<_>>();
    PotentialError::invalid_geometry(interaction, atoms, kind)
}
