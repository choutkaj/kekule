use super::*;

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Conformer {
    pub(crate) positions: Vec<Option<Point3>>,
}

impl Conformer {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_atom_capacity(atom_capacity: usize) -> Self {
        Self {
            positions: vec![None; atom_capacity],
        }
    }

    pub fn set_position(&mut self, atom: AtomId, point: Point3) {
        if self.positions.len() <= atom.index() {
            self.positions.resize(atom.index() + 1, None);
        }
        self.positions[atom.index()] = Some(point);
    }

    pub fn clear_position(&mut self, atom: AtomId) {
        if let Some(position) = self.positions.get_mut(atom.index()) {
            *position = None;
        }
    }

    pub fn position(&self, atom: AtomId) -> Option<Point3> {
        self.positions.get(atom.index()).copied().flatten()
    }

    pub fn positions(&self) -> impl Iterator<Item = (AtomId, Point3)> + '_ {
        self.positions
            .iter()
            .enumerate()
            .filter_map(|(index, point)| point.map(|point| (AtomId::new(index as u32), point)))
    }
}
