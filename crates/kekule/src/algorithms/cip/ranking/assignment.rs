use super::*;

pub(super) enum CipElementAssignment {
    Assigned(StereoDescriptor),
    Skipped(CipSkippedReason),
    Unresolved,
    Issue(CipAssignmentIssue),
}

pub(super) fn assign_cip_element(
    mol: &Molecule,
    id: StereoElementId,
    element: &StereoElement,
    options: CipAssignmentOptions,
) -> CipElementAssignment {
    if !element.is_specified() {
        return CipElementAssignment::Skipped(CipSkippedReason::UnknownConfiguration);
    }
    let assignment = match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            assign_tetrahedral_descriptor(mol, id, stereo, options)
        }
        StereoElementKind::DoubleBond(stereo) => {
            assign_double_bond_descriptor(mol, id, stereo, options)
        }
        StereoElementKind::Axis(stereo) => assign_axis_descriptor(mol, id, stereo, options),
    };
    match assignment {
        Ok(descriptor) => CipElementAssignment::Assigned(descriptor),
        Err(CipAssignmentIssue::UnresolvedPriority { .. }) => CipElementAssignment::Unresolved,
        Err(issue) => CipElementAssignment::Issue(issue),
    }
}

pub(super) fn set_stereo_descriptor(
    mol: &mut Molecule,
    id: StereoElementId,
    descriptor: StereoDescriptor,
) {
    mol.install_cip_descriptor(id, descriptor);
}

fn assign_tetrahedral_descriptor(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &TetrahedralStereo,
    options: CipAssignmentOptions,
) -> CipResult<StereoDescriptor> {
    match assign_tetrahedral_descriptor_with_auxiliary(mol, element, stereo, options, false) {
        Err(CipAssignmentIssue::UnresolvedPriority { .. }) => {
            assign_tetrahedral_descriptor_with_auxiliary(mol, element, stereo, options, true)
        }
        result => result,
    }
}

fn assign_tetrahedral_descriptor_with_auxiliary(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &TetrahedralStereo,
    options: CipAssignmentOptions,
    allow_auxiliary_descriptors: bool,
) -> CipResult<StereoDescriptor> {
    let ranked = ranked_tetrahedral_carriers(
        mol,
        element,
        stereo.center,
        &stereo.carriers,
        options,
        allow_auxiliary_descriptors,
    )?;
    tetrahedral_descriptor_from_ranked(element, stereo, &ranked)
}

pub(super) fn tetrahedral_descriptor_from_ranked(
    element: StereoElementId,
    stereo: &TetrahedralStereo,
    ranked: &RankedCarriers,
) -> CipResult<StereoDescriptor> {
    let orientation = stereo
        .orientation
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let mut priority_positions = Vec::new();
    for carrier in &ranked.carriers {
        let Some(position) = stereo
            .carriers
            .iter()
            .position(|candidate| candidate == carrier)
        else {
            return Err(CipAssignmentIssue::UnresolvedPriority { element });
        };
        priority_positions.push(position);
    }
    let even = permutation_is_even(&priority_positions);
    let descriptor_is_r = matches!(orientation, TetrahedralOrientation::Clockwise) != even;
    let descriptor = match (descriptor_is_r, ranked.pseudo_asymmetric_ordering) {
        (true, true) => StereoDescriptor::LowerR,
        (false, true) => StereoDescriptor::LowerS,
        (true, false) => StereoDescriptor::R,
        (false, false) => StereoDescriptor::S,
    };
    Ok(descriptor)
}

fn assign_double_bond_descriptor(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &DoubleBondStereo,
    options: CipAssignmentOptions,
) -> CipResult<StereoDescriptor> {
    let orientation = stereo
        .orientation
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let left_carriers = double_bond_endpoint_carriers(mol, stereo.left, stereo.right, stereo.bond);
    let right_carriers = double_bond_endpoint_carriers(mol, stereo.right, stereo.left, stereo.bond);
    let left_ranked = ranked_carriers(
        mol,
        element,
        stereo.left,
        &left_carriers,
        options,
        true,
        false,
    )?;
    let right_ranked = ranked_carriers(
        mol,
        element,
        stereo.right,
        &right_carriers,
        options,
        true,
        false,
    )?;
    double_bond_descriptor_from_ranked(
        element,
        orientation,
        (stereo.left_carrier, stereo.right_carrier),
        &left_ranked,
        &right_ranked,
    )
}

fn assign_axis_descriptor(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &AxisStereo,
    options: CipAssignmentOptions,
) -> CipResult<StereoDescriptor> {
    let orientation = stereo
        .orientation
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let bond = mol
        .bond(stereo.axis)
        .map_err(|_| CipAssignmentIssue::UnresolvedPriority { element })?;
    let (left, right) = bond.endpoints();
    let (left_reference, right_reference) =
        axis_reference_carriers(mol, element, stereo, left, right)?;
    let left_ranked = ranked_carriers(
        mol,
        element,
        left,
        &axis_endpoint_carriers(mol, left, right, stereo.axis),
        options,
        true,
        true,
    )?;
    let right_ranked = ranked_carriers(
        mol,
        element,
        right,
        &axis_endpoint_carriers(mol, right, left, stereo.axis),
        options,
        true,
        true,
    )?;
    axis_descriptor_from_ranked(
        element,
        orientation,
        (left_reference, right_reference),
        &left_ranked,
        &right_ranked,
    )
}

fn ranked_endpoint_relation(
    element: StereoElementId,
    (left_reference, right_reference): (StereoCarrier, StereoCarrier),
    left: &RankedCarriers,
    right: &RankedCarriers,
) -> CipResult<(bool, bool)> {
    let left_top = left
        .carriers
        .first()
        .copied()
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let right_top = right
        .carriers
        .first()
        .copied()
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    Ok((
        (left_reference != left_top) != (right_reference != right_top),
        left.pseudo_asymmetric_ordering != right.pseudo_asymmetric_ordering,
    ))
}

pub(super) fn double_bond_descriptor_from_ranked(
    element: StereoElementId,
    mut orientation: DoubleBondOrientation,
    references: (StereoCarrier, StereoCarrier),
    left: &RankedCarriers,
    right: &RankedCarriers,
) -> CipResult<StereoDescriptor> {
    let (inverted, pseudo) = ranked_endpoint_relation(element, references, left, right)?;
    if inverted {
        orientation = orientation.inverted();
    }
    Ok(match (orientation, pseudo) {
        (DoubleBondOrientation::Together, true) => StereoDescriptor::SeqCis,
        (DoubleBondOrientation::Opposite, true) => StereoDescriptor::SeqTrans,
        (DoubleBondOrientation::Together, false) => StereoDescriptor::Z,
        (DoubleBondOrientation::Opposite, false) => StereoDescriptor::E,
    })
}

pub(super) fn axis_descriptor_from_ranked(
    element: StereoElementId,
    mut orientation: AxisOrientation,
    references: (StereoCarrier, StereoCarrier),
    left: &RankedCarriers,
    right: &RankedCarriers,
) -> CipResult<StereoDescriptor> {
    let (inverted, pseudo) = ranked_endpoint_relation(element, references, left, right)?;
    if inverted {
        orientation = orientation.inverted();
    }
    Ok(match (orientation, pseudo) {
        (AxisOrientation::CounterClockwise, true) => StereoDescriptor::LowerM,
        (AxisOrientation::Clockwise, true) => StereoDescriptor::LowerP,
        (AxisOrientation::CounterClockwise, false) => StereoDescriptor::M,
        (AxisOrientation::Clockwise, false) => StereoDescriptor::P,
    })
}

pub(super) fn axis_reference_carriers(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &AxisStereo,
    left: AtomId,
    right: AtomId,
) -> CipResult<(StereoCarrier, StereoCarrier)> {
    if stereo.carriers.len() != 2 {
        return Err(CipAssignmentIssue::UnresolvedPriority { element });
    }
    let mut left_reference = None;
    let mut right_reference = None;
    for carrier in &stereo.carriers {
        let StereoCarrier::Atom(atom) = carrier else {
            return Err(CipAssignmentIssue::UnresolvedPriority { element });
        };
        let adjacent_left = mol.bond_between(left, *atom).ok().flatten().is_some();
        let adjacent_right = mol.bond_between(right, *atom).ok().flatten().is_some();
        match (adjacent_left, adjacent_right) {
            (true, false) if left_reference.is_none() => left_reference = Some(*carrier),
            (false, true) if right_reference.is_none() => right_reference = Some(*carrier),
            _ => return Err(CipAssignmentIssue::UnresolvedPriority { element }),
        }
    }
    match (left_reference, right_reference) {
        (Some(left), Some(right)) => Ok((left, right)),
        _ => Err(CipAssignmentIssue::UnresolvedPriority { element }),
    }
}

pub(super) fn axis_endpoint_carriers(
    mol: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    axis: BondId,
) -> Vec<StereoCarrier> {
    let mut carriers = Vec::new();
    if let Ok(incident) = mol.incident_bonds(endpoint) {
        for (bond_id, bond) in incident {
            if bond_id != axis {
                carriers.push(StereoCarrier::Atom(bond.other_atom(endpoint)));
            }
        }
    }
    if atom_hydrogen_count(mol, endpoint) > 0
        && mol
            .bond_between(endpoint, other_endpoint)
            .ok()
            .flatten()
            .is_some()
    {
        carriers.push(StereoCarrier::ImplicitHydrogen);
    }
    carriers
}

fn ranked_carriers(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    carriers: &[StereoCarrier],
    options: CipAssignmentOptions,
    allow_auxiliary_descriptors: bool,
    atropisomer_mode: bool,
) -> CipResult<RankedCarriers> {
    if allow_auxiliary_descriptors {
        let constitutional = carrier_signatures(
            mol,
            element,
            root,
            carriers,
            options,
            false,
            atropisomer_mode,
        )?;
        if let Ok(ranked) = rank_carrier_signatures(element, &constitutional, None) {
            return Ok(ranked);
        }
    }
    let signatures = carrier_signatures(
        mol,
        element,
        root,
        carriers,
        options,
        allow_auxiliary_descriptors,
        atropisomer_mode,
    )?;
    rank_carrier_signatures(element, &signatures, None)
}

fn ranked_tetrahedral_carriers(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    carriers: &[StereoCarrier],
    options: CipAssignmentOptions,
    allow_auxiliary_descriptors: bool,
) -> CipResult<RankedCarriers> {
    let signatures = carrier_signatures(
        mol,
        element,
        root,
        carriers,
        options,
        allow_auxiliary_descriptors,
        false,
    )?;
    match rank_carrier_signatures(element, &signatures, None) {
        Ok(ranked) => Ok(ranked),
        Err(CipAssignmentIssue::UnresolvedPriority { .. })
            if allow_auxiliary_descriptors && carriers.len() == 4 =>
        {
            rank_tetrahedral_signatures_with_rule6(element, &signatures)
        }
        Err(issue) => Err(issue),
    }
}

pub(super) fn element_is_finally_nonstereogenic(
    mol: &Molecule,
    element: StereoElementId,
    stereo_element: &StereoElement,
    options: CipAssignmentOptions,
) -> CipResult<bool> {
    if !stereo_element.is_specified() {
        return Ok(false);
    }
    match &stereo_element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            tetrahedral_final_tie_is_nonstereogenic(mol, element, stereo, options)
        }
        StereoElementKind::DoubleBond(stereo) => {
            double_bond_final_tie_is_nonstereogenic(mol, element, stereo, options)
        }
        StereoElementKind::Axis(stereo) => {
            axis_final_tie_is_nonstereogenic(mol, element, stereo, options)
        }
    }
}

fn tetrahedral_final_tie_is_nonstereogenic(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &TetrahedralStereo,
    options: CipAssignmentOptions,
) -> CipResult<bool> {
    let signatures = carrier_signatures(
        mol,
        element,
        stereo.center,
        &stereo.carriers,
        options,
        true,
        false,
    )?;
    match rank_tetrahedral_signatures_with_rule6(element, &signatures) {
        Ok(_) => Ok(false),
        Err(CipAssignmentIssue::UnresolvedPriority { .. }) => {
            Ok(grouped_signature_indices(&signatures)
                .iter()
                .any(|group| group.len() > 1))
        }
        Err(issue) => Err(issue),
    }
}

fn double_bond_final_tie_is_nonstereogenic(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &DoubleBondStereo,
    options: CipAssignmentOptions,
) -> CipResult<bool> {
    let left_carriers = double_bond_endpoint_carriers(mol, stereo.left, stereo.right, stereo.bond);
    let right_carriers = double_bond_endpoint_carriers(mol, stereo.right, stereo.left, stereo.bond);
    Ok(endpoint_final_tie_is_nonstereogenic(
        mol,
        element,
        stereo.left,
        &left_carriers,
        options,
        false,
    )? || endpoint_final_tie_is_nonstereogenic(
        mol,
        element,
        stereo.right,
        &right_carriers,
        options,
        false,
    )?)
}

fn axis_final_tie_is_nonstereogenic(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &AxisStereo,
    options: CipAssignmentOptions,
) -> CipResult<bool> {
    let bond = mol
        .bond(stereo.axis)
        .map_err(|_| CipAssignmentIssue::UnresolvedPriority { element })?;
    let (left, right) = bond.endpoints();
    axis_reference_carriers(mol, element, stereo, left, right)?;
    let left_carriers = axis_endpoint_carriers(mol, left, right, stereo.axis);
    let right_carriers = axis_endpoint_carriers(mol, right, left, stereo.axis);
    Ok(
        endpoint_final_tie_is_nonstereogenic(mol, element, left, &left_carriers, options, true)?
            || endpoint_final_tie_is_nonstereogenic(
                mol,
                element,
                right,
                &right_carriers,
                options,
                true,
            )?,
    )
}

fn endpoint_final_tie_is_nonstereogenic(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    carriers: &[StereoCarrier],
    options: CipAssignmentOptions,
    atropisomer_mode: bool,
) -> CipResult<bool> {
    if carriers.len() < 2 {
        return Ok(false);
    }
    let signatures = carrier_signatures(
        mol,
        element,
        root,
        carriers,
        options,
        true,
        atropisomer_mode,
    )?;
    match rank_carrier_signatures(element, &signatures, None) {
        Ok(_) => Ok(false),
        Err(CipAssignmentIssue::UnresolvedPriority { .. }) => {
            Ok(grouped_signature_indices(&signatures)
                .iter()
                .any(|group| group.len() > 1))
        }
        Err(issue) => Err(issue),
    }
}

fn carrier_signatures(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    carriers: &[StereoCarrier],
    options: CipAssignmentOptions,
    allow_auxiliary_descriptors: bool,
    atropisomer_mode: bool,
) -> CipResult<Vec<(StereoCarrier, LigandSignature)>> {
    let atomic_number_fractions = cip_atomic_number_fractions(mol);
    if allow_auxiliary_descriptors
        && mol
            .stereo_elements()
            .any(|(id, stereo)| id != element && stereo.is_specified())
    {
        let mut descriptor_context = DescriptorContext::new(element);
        let aux_graph = build_auxiliary_graph(
            mol,
            element,
            root,
            options,
            &atomic_number_fractions,
            atropisomer_mode,
        )?;
        precompute_auxiliary_descriptors(
            mol,
            &mut descriptor_context,
            &aux_graph,
            options,
            &atomic_number_fractions,
            atropisomer_mode,
        );
        let build_context = LigandBuildContext {
            mol,
            element,
            descriptor_context: &descriptor_context,
            options,
            atomic_number_fractions: &atomic_number_fractions,
            atropisomer_mode,
        };
        let signatures = build_carrier_signatures(&build_context, root, carriers)?;
        return Ok(signatures);
    }
    let descriptor_context = DescriptorContext::new(element);
    let mut build_context = LigandBuildContext {
        mol,
        element,
        descriptor_context: &descriptor_context,
        options,
        atomic_number_fractions: &atomic_number_fractions,
        atropisomer_mode,
    };
    let mut depth = 0;
    loop {
        build_context.options.max_depth = depth;
        let signatures = build_carrier_signatures(&build_context, root, carriers)?;
        if signatures.iter().all(|(_, signature)| !signature.truncated)
            || rank_carrier_signatures(element, &signatures, None).is_ok()
        {
            return Ok(signatures);
        }
        if depth == options.max_depth {
            return Err(CipAssignmentIssue::DepthLimitExceeded {
                element,
                max_depth: options.max_depth,
            });
        }
        depth = depth.saturating_mul(2).max(1).min(options.max_depth);
    }
}

fn build_carrier_signatures(
    context: &LigandBuildContext<'_>,
    root: AtomId,
    carriers: &[StereoCarrier],
) -> CipResult<Vec<(StereoCarrier, LigandSignature)>> {
    carriers
        .iter()
        .copied()
        .map(|carrier| {
            carrier_signature(context, carrier, root).map(|signature| (carrier, signature))
        })
        .collect::<CipResult<Vec<_>>>()
}

pub(super) fn rank_carrier_signatures(
    element: StereoElementId,
    signatures: &[(StereoCarrier, LigandSignature)],
    rule6_reference: Option<AtomId>,
) -> CipResult<RankedCarriers> {
    let mut pseudo_asymmetric_pair_count = 0usize;
    for left in 0..signatures.len() {
        for right in (left + 1)..signatures.len() {
            let comparison = signatures[left]
                .1
                .compare_with_rule6_reference(&signatures[right].1, rule6_reference);
            if comparison.ordering == Ordering::Equal {
                return Err(CipAssignmentIssue::UnresolvedPriority { element });
            }
            if comparison.pseudo_asymmetric {
                pseudo_asymmetric_pair_count += 1;
            }
        }
    }
    let mut signatures = signatures.to_vec();
    signatures.sort_by(|left, right| {
        right
            .1
            .compare_with_rule6_reference(&left.1, rule6_reference)
            .ordering
    });
    Ok(RankedCarriers {
        carriers: signatures.into_iter().map(|(carrier, _)| carrier).collect(),
        pseudo_asymmetric_ordering: pseudo_asymmetric_pair_count == 1,
    })
}

pub(super) fn rank_tetrahedral_signatures_with_rule6(
    element: StereoElementId,
    signatures: &[(StereoCarrier, LigandSignature)],
) -> CipResult<RankedCarriers> {
    let groups = grouped_signature_indices(signatures);
    match groups.len() {
        2 => {
            let Some(reference_index) = groups.iter().flatten().copied().nth(1) else {
                return Err(CipAssignmentIssue::UnresolvedPriority { element });
            };
            let Some(reference) = carrier_rule6_atom(signatures[reference_index].0) else {
                return Err(CipAssignmentIssue::UnresolvedPriority { element });
            };
            let ranked = rank_carrier_signatures(element, signatures, Some(reference))?;
            reject_rule6_parity_unstable_references(element, signatures, &ranked)?;
            Ok(ranked)
        }
        1 => rank_s4_tetrahedral_signatures_with_rule6(element, signatures, &groups[0]),
        _ => Err(CipAssignmentIssue::UnresolvedPriority { element }),
    }
}

fn rank_s4_tetrahedral_signatures_with_rule6(
    element: StereoElementId,
    signatures: &[(StereoCarrier, LigandSignature)],
    group: &[usize],
) -> CipResult<RankedCarriers> {
    let mut stable_ranking: Option<RankedCarriers> = None;
    for index in group {
        let Some(reference) = carrier_rule6_atom(signatures[*index].0) else {
            continue;
        };
        let ranking = match rank_carrier_signatures(element, signatures, Some(reference)) {
            Ok(ranking) => ranking,
            Err(CipAssignmentIssue::UnresolvedPriority { .. }) => continue,
            Err(issue) => return Err(issue),
        };
        if let Some(stable) = &stable_ranking {
            if carrier_permutation_is_odd(&stable.carriers, &ranking.carriers).unwrap_or(true) {
                return Err(CipAssignmentIssue::UnresolvedPriority { element });
            }
        } else {
            stable_ranking = Some(ranking);
        }
    }
    stable_ranking.ok_or(CipAssignmentIssue::UnresolvedPriority { element })
}

fn reject_rule6_parity_unstable_references(
    element: StereoElementId,
    signatures: &[(StereoCarrier, LigandSignature)],
    stable: &RankedCarriers,
) -> CipResult<()> {
    for (carrier, _) in signatures {
        let Some(reference) = carrier_rule6_atom(*carrier) else {
            continue;
        };
        let ranking = match rank_carrier_signatures(element, signatures, Some(reference)) {
            Ok(ranking) => ranking,
            Err(CipAssignmentIssue::UnresolvedPriority { .. }) => continue,
            Err(issue) => return Err(issue),
        };
        if carrier_permutation_is_odd(&stable.carriers, &ranking.carriers).unwrap_or(true) {
            return Err(CipAssignmentIssue::UnresolvedPriority { element });
        }
    }
    Ok(())
}

fn grouped_signature_indices(signatures: &[(StereoCarrier, LigandSignature)]) -> Vec<Vec<usize>> {
    let mut indices = (0..signatures.len()).collect::<Vec<_>>();
    indices.sort_by(|left, right| signatures[*right].1.compare(&signatures[*left].1));
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for index in indices {
        if let Some(last) = groups.last_mut() {
            if signatures[last[0]].1.compare(&signatures[index].1) == Ordering::Equal {
                last.push(index);
                continue;
            }
        }
        groups.push(vec![index]);
    }
    groups
}

fn carrier_rule6_atom(carrier: StereoCarrier) -> Option<AtomId> {
    match carrier {
        StereoCarrier::Atom(atom) => Some(atom),
        StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => None,
    }
}

fn carrier_permutation_is_odd(left: &[StereoCarrier], right: &[StereoCarrier]) -> Option<bool> {
    if left.len() != right.len() {
        return None;
    }
    let mut positions = Vec::with_capacity(left.len());
    for carrier in left {
        positions.push(right.iter().position(|candidate| candidate == carrier)?);
    }
    Some(!permutation_is_even(&positions))
}
