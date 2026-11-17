//! Joining a frame table into a structure set by step:
//! [`StructureSet::join_frames`] and [`JoinReport`].

use crate::properties::{LiveRows, PropertyColumn, PropertyKey, PropertyTableMut};
use crate::structure::{
    GeometryIndex, GeometryStorage, Requirement, StructureSet, StructureSetError,
};

use super::header::STEP;
use super::read::FrameTable;

/// What [`StructureSet::join_frames`] joined.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JoinReport {
    columns: Vec<PropertyKey>,
    matched: usize,
    unmatched: usize,
}

impl JoinReport {
    /// The keys of the joined columns, in the table's order.
    pub fn columns(&self) -> &[PropertyKey] {
        &self.columns
    }

    /// The geometries whose step has a row in the table.
    pub const fn matched(&self) -> usize {
        self.matched
    }

    /// The geometries whose step has no row, or that have no step; their
    /// joined cells are missing.
    pub const fn unmatched(&self) -> usize {
        self.unmatched
    }
}

impl<S: GeometryStorage> StructureSet<S> {
    /// Joins the real and integer columns of a frame table into the set's
    /// columns by step, and installs the table's sampling record.
    ///
    /// Each geometry takes the row of its step: the step of its order
    /// ([`Self::steps`]), or else of its integer column `step`. A geometry
    /// whose step has no row, or that has no step, gets missing cells, and a
    /// row without a geometry is left out. The table's step and time columns
    /// and its boolean and text columns are not joined. A joined column
    /// replaces the set's column of the same key, so joining a table again
    /// after it grew fills the cells its new rows hold.
    ///
    /// The set takes the table's sampling record; it keeps its own when the
    /// table has none. Failure changes nothing.
    ///
    /// ```
    /// use kekule::properties::PropertyKey;
    /// use kekule::structure::{Geometry, Order, Positions, StructureSet};
    /// use kekule::traj::frames::read_frame_table_from;
    ///
    /// let topology = kekule::smiles::to_topology("O")?;
    /// let geometry = || Geometry::new(Positions::zeros(topology.atom_count()));
    /// let mut set = StructureSet::from_geometries(topology.clone(), [geometry(), geometry()])?;
    /// set.set_order(Some(Order::Step(vec![0, 1000])))?;
    /// let text = "#kekule-frames\t1\n#column\tstep\tint\n#column\tbias\treal\tkJ/mol\n\
    ///             step\tbias\n0\t1.5\n500\t2\n1000\t2.5\n";
    /// let report = set.join_frames(&read_frame_table_from(text.as_bytes())?)?;
    /// assert_eq!(report.matched(), 2);
    /// let bias = set.columns().get(&PropertyKey::new("bias")?).unwrap();
    /// assert_eq!(bias.value(1)?, Some(kekule::properties::PropertyValue::Real {
    ///     value: 2.5,
    ///     unit: kekule::units::KILOJOULE_PER_MOLE,
    /// }));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// `Missing(Steps)` for a set without steps or an integer column
    /// `step`, [`StructureSetError::SamplingMismatch`] when the set has
    /// another sampling record than the table, and
    /// [`StructureSetError::SamplingColumn`] when the set's own record would
    /// name a column that holds other values than energies.
    pub fn join_frames(&mut self, table: &FrameTable) -> Result<JoinReport, StructureSetError> {
        let sampling = match (self.sampling(), table.header().sampling()) {
            (Some(own), Some(table)) if own != table => {
                return Err(StructureSetError::SamplingMismatch)
            }
            (own, table) => table.or(own).cloned(),
        };
        let rows = self.joined_rows(table)?;
        let mut columns = self.columns().clone();
        let mut joined = Vec::new();
        let mut guard = PropertyTableMut::<GeometryIndex>::new(&mut columns, LiveRows::All, None);
        for spec in table.header().columns() {
            let key = spec.key();
            if Some(key) == table.header().time() {
                continue;
            }
            let column = match table.column(key) {
                Some(PropertyColumn::Real { unit, values }) => PropertyColumn::Real {
                    unit: *unit,
                    values: rows
                        .iter()
                        .map(|row| row.and_then(|row| values[row]))
                        .collect(),
                },
                Some(PropertyColumn::Int(values)) => PropertyColumn::Int(
                    rows.iter()
                        .map(|row| row.and_then(|row| values[row]))
                        .collect(),
                ),
                _ => continue,
            };
            guard.remove(key);
            guard.insert(key.clone(), column)?;
            joined.push(key.clone());
        }
        if let Some(sampling) = &sampling {
            sampling.check_columns(&columns)?;
        }
        let matched = rows.iter().flatten().count();
        self.set_columns(columns);
        self.set_sampling(sampling)
            .expect("the record's columns are checked");
        Ok(JoinReport {
            columns: joined,
            matched,
            unmatched: rows.len() - matched,
        })
    }

    /// The table row of each geometry's step, when it has one.
    fn joined_rows(&self, table: &FrameTable) -> Result<Vec<Option<usize>>, StructureSetError> {
        if let Some(steps) = self.steps() {
            return Ok(steps.iter().map(|&step| table.row_of_step(step)).collect());
        }
        let key = PropertyKey::new(STEP).expect("\"step\" is a valid key");
        match self.columns().get(&key) {
            Some(PropertyColumn::Int(steps)) => Ok(steps
                .iter()
                .map(|step| {
                    step.and_then(|step| u64::try_from(step).ok())
                        .and_then(|step| table.row_of_step(step))
                })
                .collect()),
            _ => Err(Requirement::Steps.into()),
        }
    }
}
