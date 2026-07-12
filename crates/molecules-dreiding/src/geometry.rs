use molecules::core::Point3;

const MIN_NORM_SQ: f64 = 1.0e-24;

#[derive(Debug, Clone, Copy)]
pub(crate) enum GeometryError {
    Coincident,
    DegenerateAngle,
    DegenerateDihedral,
    DegenerateInversion,
}

#[derive(Clone, Copy)]
struct Dual<const N: usize> {
    value: f64,
    derivative: [f64; N],
}

impl<const N: usize> Dual<N> {
    fn variable(value: f64, index: usize) -> Self {
        let mut derivative = [0.0; N];
        derivative[index] = 1.0;
        Self { value, derivative }
    }

    fn add(self, other: Self) -> Self {
        let mut derivative = [0.0; N];
        for (index, value) in derivative.iter_mut().enumerate() {
            *value = self.derivative[index] + other.derivative[index];
        }
        Self {
            value: self.value + other.value,
            derivative,
        }
    }

    fn sub(self, other: Self) -> Self {
        self.add(other.scale(-1.0))
    }

    fn mul(self, other: Self) -> Self {
        let mut derivative = [0.0; N];
        for (index, value) in derivative.iter_mut().enumerate() {
            *value = self.derivative[index] * other.value + self.value * other.derivative[index];
        }
        Self {
            value: self.value * other.value,
            derivative,
        }
    }

    fn scale(self, scale: f64) -> Self {
        let mut derivative = self.derivative;
        for value in &mut derivative {
            *value *= scale;
        }
        Self {
            value: self.value * scale,
            derivative,
        }
    }

    fn reciprocal(self) -> Self {
        let inverse = self.value.recip();
        let mut derivative = self.derivative;
        let scale = -inverse * inverse;
        for value in &mut derivative {
            *value *= scale;
        }
        Self {
            value: inverse,
            derivative,
        }
    }

    fn div(self, other: Self) -> Self {
        self.mul(other.reciprocal())
    }

    fn sqrt(self) -> Self {
        let root = self.value.sqrt();
        let mut derivative = self.derivative;
        let scale = 0.5 / root;
        for value in &mut derivative {
            *value *= scale;
        }
        Self {
            value: root,
            derivative,
        }
    }
}

#[derive(Clone, Copy)]
struct DualVector<const N: usize> {
    x: Dual<N>,
    y: Dual<N>,
    z: Dual<N>,
}

impl<const N: usize> DualVector<N> {
    fn from_point(point: Point3, slot: usize) -> Self {
        Self {
            x: Dual::variable(point.x, slot * 3),
            y: Dual::variable(point.y, slot * 3 + 1),
            z: Dual::variable(point.z, slot * 3 + 2),
        }
    }

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x.sub(other.x),
            y: self.y.sub(other.y),
            z: self.z.sub(other.z),
        }
    }

    fn scale(self, scale: Dual<N>) -> Self {
        Self {
            x: self.x.mul(scale),
            y: self.y.mul(scale),
            z: self.z.mul(scale),
        }
    }

    fn dot(self, other: Self) -> Dual<N> {
        self.x
            .mul(other.x)
            .add(self.y.mul(other.y))
            .add(self.z.mul(other.z))
    }

    fn cross(self, other: Self) -> Self {
        Self {
            x: self.y.mul(other.z).sub(self.z.mul(other.y)),
            y: self.z.mul(other.x).sub(self.x.mul(other.z)),
            z: self.x.mul(other.y).sub(self.y.mul(other.x)),
        }
    }

    fn norm_sq(self) -> Dual<N> {
        self.dot(self)
    }

    fn normalized(self, error: GeometryError) -> Result<Self, GeometryError> {
        let norm_sq = self.norm_sq();
        if norm_sq.value <= MIN_NORM_SQ {
            return Err(error);
        }
        Ok(self.scale(norm_sq.sqrt().reciprocal()))
    }
}

pub(crate) fn displacement(
    first: Point3,
    second: Point3,
) -> Result<([f64; 3], f64), GeometryError> {
    let vector = [first.x - second.x, first.y - second.y, first.z - second.z];
    let norm_sq = vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
    if norm_sq <= MIN_NORM_SQ {
        return Err(GeometryError::Coincident);
    }
    Ok((vector, norm_sq))
}

pub(crate) fn angle_cosine(
    first: Point3,
    center: Point3,
    third: Point3,
) -> Result<(f64, [f64; 9]), GeometryError> {
    let first = DualVector::<9>::from_point(first, 0);
    let center = DualVector::<9>::from_point(center, 1);
    let third = DualVector::<9>::from_point(third, 2);
    let left = first
        .sub(center)
        .normalized(GeometryError::DegenerateAngle)?;
    let right = third
        .sub(center)
        .normalized(GeometryError::DegenerateAngle)?;
    let cosine = left.dot(right);
    Ok((cosine.value.clamp(-1.0, 1.0), cosine.derivative))
}

pub(crate) fn torsion(points: [Point3; 4]) -> Result<(f64, f64, [f64; 12]), GeometryError> {
    let p0 = DualVector::<12>::from_point(points[0], 0);
    let p1 = DualVector::<12>::from_point(points[1], 1);
    let p2 = DualVector::<12>::from_point(points[2], 2);
    let p3 = DualVector::<12>::from_point(points[3], 3);

    let b0 = p0.sub(p1);
    let b1 = p2.sub(p1);
    let b2 = p3.sub(p2);
    let b1_norm_sq = b1.norm_sq();
    if b1_norm_sq.value <= MIN_NORM_SQ {
        return Err(GeometryError::DegenerateDihedral);
    }
    let b1_unit = b1.normalized(GeometryError::DegenerateDihedral)?;
    let v = b0.sub(b1.scale(b0.dot(b1).div(b1_norm_sq)));
    let w = b2.sub(b1.scale(b2.dot(b1).div(b1_norm_sq)));
    let v_unit = v.normalized(GeometryError::DegenerateDihedral)?;
    let w_unit = w.normalized(GeometryError::DegenerateDihedral)?;
    let cosine = v_unit.dot(w_unit);
    let sine = b1_unit.cross(v_unit).dot(w_unit);
    let mut derivative = [0.0; 12];
    for (index, value) in derivative.iter_mut().enumerate() {
        *value = cosine.value * sine.derivative[index] - sine.value * cosine.derivative[index];
    }
    Ok((
        cosine.value.clamp(-1.0, 1.0),
        sine.value.clamp(-1.0, 1.0),
        derivative,
    ))
}

pub(crate) fn inversion_cosine(points: [Point3; 4]) -> Result<(f64, [f64; 12]), GeometryError> {
    let center = DualVector::<12>::from_point(points[0], 0);
    let axis = DualVector::<12>::from_point(points[1], 1);
    let plane1 = DualVector::<12>::from_point(points[2], 2);
    let plane2 = DualVector::<12>::from_point(points[3], 3);
    let normal = plane1
        .sub(center)
        .cross(plane2.sub(center))
        .normalized(GeometryError::DegenerateInversion)?;
    let axis = axis
        .sub(center)
        .normalized(GeometryError::DegenerateInversion)?;
    let cosine = normal.dot(axis);
    Ok((cosine.value.clamp(-1.0, 1.0), cosine.derivative))
}

pub(crate) fn hydrogen_bond_cosine(
    donor: Point3,
    hydrogen: Point3,
    acceptor: Point3,
) -> Result<(f64, [f64; 9]), GeometryError> {
    let donor = DualVector::<9>::from_point(donor, 0);
    let hydrogen = DualVector::<9>::from_point(hydrogen, 1);
    let acceptor = DualVector::<9>::from_point(acceptor, 2);
    let donor_leg = donor
        .sub(hydrogen)
        .normalized(GeometryError::DegenerateAngle)?;
    let acceptor_leg = acceptor
        .sub(hydrogen)
        .normalized(GeometryError::DegenerateAngle)?;
    // dreid-kernel expects 1.0 for a linear D-H...A arrangement.
    let cosine = donor_leg.dot(acceptor_leg).scale(-1.0);
    Ok((cosine.value.clamp(-1.0, 1.0), cosine.derivative))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_hydrogen_bond_has_positive_cosine() {
        let (cosine, _) = hydrogen_bond_cosine(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        )
        .unwrap();
        assert!((cosine - 1.0).abs() < 1.0e-12);
    }
}
