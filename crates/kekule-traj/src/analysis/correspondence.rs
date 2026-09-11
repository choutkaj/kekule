//! Independent-reference workflows, sharing the existing frame and RMSD kernels.
use super::*;

impl Trajectory {
    /// Returns a copy fitted to an independent model through explicit atom pairs.
    ///
    /// Only the paired atoms determine the fit; all moving positions, vectors,
    /// and cells are transformed. The original moving topology and annotations
    /// are retained, even when the reference has a different atom count or order.
    pub fn superpose_to_model(
        &self,
        reference: ModelView<'_>,
        correspondence: &AtomCorrespondence,
    ) -> Result<Self, SuperpositionError> {
        self.superpose_to_model_with_options(reference, correspondence, KabschOptions::default())
    }

    /// Fits to an independent reference with weights in correspondence pair order.
    pub fn superpose_to_model_with_options(
        &self,
        reference: ModelView<'_>,
        correspondence: &AtomCorrespondence,
        options: KabschOptions<'_>,
    ) -> Result<Self, SuperpositionError> {
        let frames = self.model_superposition_frames(reference, correspondence, options)?;
        self.with_frames(frames)
            .map_err(|source| SuperpositionError::TrajectoryPublication(Box::new(source)))
    }

    /// Fits to an independent reference transactionally, leaving this trajectory intact on failure.
    pub fn superpose_to_model_in_place(
        &mut self,
        reference: ModelView<'_>,
        correspondence: &AtomCorrespondence,
    ) -> Result<(), SuperpositionError> {
        self.superpose_to_model_in_place_with_options(
            reference,
            correspondence,
            KabschOptions::default(),
        )
    }

    /// Transactional fitting with explicit weights and periodic policy.
    pub fn superpose_to_model_in_place_with_options(
        &mut self,
        reference: ModelView<'_>,
        correspondence: &AtomCorrespondence,
        options: KabschOptions<'_>,
    ) -> Result<(), SuperpositionError> {
        let frames = self.model_superposition_frames(reference, correspondence, options)?;
        self.replace_frames(frames)
            .map_err(|source| SuperpositionError::TrajectoryPublication(Box::new(source)))
    }

    fn model_superposition_frames(
        &self,
        reference: ModelView<'_>,
        correspondence: &AtomCorrespondence,
        options: KabschOptions<'_>,
    ) -> Result<Vec<TrajectoryFrame>, SuperpositionError> {
        correspondence
            .ensure_compatible(self.topology(), reference.topology())
            .map_err(SuperpositionError::Correspondence)?;
        let superposer =
            FrameSuperposer::with_correspondence_and_options(reference, correspondence, options);
        self.frames()
            .enumerate()
            .map(|(index, frame)| superposer.superpose(index, frame))
            .collect()
    }

    /// Measures direct RMSD to an independent model without fitting or mutation.
    ///
    /// Only the explicitly paired atoms are measured, in pair order. The two
    /// systems may differ outside these atoms. Results use canonical length units.
    pub fn rmsd_to_model(
        &self,
        reference: ModelView<'_>,
        correspondence: &AtomCorrespondence,
    ) -> Result<Quantity<Vec<f64>>, RmsdError> {
        self.rmsd_to_model_with_options(reference, correspondence, RmsdOptions::default())
    }

    /// Measures paired atoms with explicit weights and periodic policy.
    pub fn rmsd_to_model_with_options(
        &self,
        reference: ModelView<'_>,
        correspondence: &AtomCorrespondence,
        options: RmsdOptions<'_>,
    ) -> Result<Quantity<Vec<f64>>, RmsdError> {
        self.validate_measurement_correspondence(reference, correspondence)?;
        let weights = NormalizedRmsdWeights::new(options.weighting, correspondence.len())?;
        let mut values = Vec::with_capacity(self.len());
        for (frame, moving) in self.frames().enumerate() {
            if options.periodic_policy == PeriodicRmsdPolicy::RejectPeriodic
                && (moving.cell().is_some() || reference.cell().is_some())
            {
                return Err(RmsdError::PeriodicCoordinates {
                    frame,
                    moving: moving.cell().is_some(),
                    reference: reference.cell().is_some(),
                });
            }
            values.push(measure_rmsd(
                moving.as_model(),
                reference,
                correspondence.index_pairs().iter().copied(),
                weights,
                None,
                frame,
            )?);
        }
        Ok(Quantity::new(values, CANONICAL_LENGTH_UNIT))
    }

    /// Fits one set of pairs and measures a separate set against an independent model.
    ///
    /// For example, pair backbone atoms for fitting and ligand atoms for measuring
    /// motion. No transformed frame copies are materialized. Both correspondences
    /// must bind this trajectory's topology to the reference topology.
    pub fn aligned_rmsd_to_model(
        &self,
        reference: ModelView<'_>,
        fit: &AtomCorrespondence,
        measurement: &AtomCorrespondence,
    ) -> Result<Quantity<Vec<f64>>, RmsdError> {
        self.aligned_rmsd_to_model_with_options(
            reference,
            fit,
            measurement,
            AlignedRmsdOptions::default(),
        )
    }

    /// Fits and measures explicit pairs with separate fit and measurement weights.
    pub fn aligned_rmsd_to_model_with_options(
        &self,
        reference: ModelView<'_>,
        fit: &AtomCorrespondence,
        measurement: &AtomCorrespondence,
        options: AlignedRmsdOptions<'_>,
    ) -> Result<Quantity<Vec<f64>>, RmsdError> {
        fit.ensure_compatible(self.topology(), reference.topology())
            .map_err(RmsdError::Correspondence)?;
        self.validate_measurement_correspondence(reference, measurement)?;
        let weights = NormalizedRmsdWeights::new(options.measurement_weighting, measurement.len())?;
        let mut values = Vec::with_capacity(self.len());
        for (frame, moving) in self.frames().enumerate() {
            let alignment = kabsch_with_correspondence_and_options(
                moving.as_model(),
                reference,
                fit,
                options.superposition,
            )
            .map_err(|source| RmsdError::Alignment { frame, source })?;
            values.push(measure_rmsd(
                moving.as_model(),
                reference,
                measurement.index_pairs().iter().copied(),
                weights,
                Some(alignment.transform()),
                frame,
            )?);
        }
        Ok(Quantity::new(values, CANONICAL_LENGTH_UNIT))
    }

    fn validate_measurement_correspondence(
        &self,
        reference: ModelView<'_>,
        correspondence: &AtomCorrespondence,
    ) -> Result<(), RmsdError> {
        correspondence
            .ensure_compatible(self.topology(), reference.topology())
            .map_err(RmsdError::Correspondence)?;
        if correspondence.is_empty() {
            return Err(RmsdError::EmptySelection);
        }
        Ok(())
    }
}
