use super::*;

pub(super) fn descriptor_is_absolute_tetrahedral(descriptor: StereoDescriptor) -> bool {
    matches!(descriptor, StereoDescriptor::R | StereoDescriptor::S)
}

pub(super) enum CipElementAssignment {
    Assigned(StereoDescriptor),
    Skipped(CipSkippedReason),
    Deferred,
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
        StereoElementKind::DoubleBond(stereo) => match double_bond_cip_stereogenic(mol, stereo) {
            Some(false) => return CipElementAssignment::Skipped(CipSkippedReason::NotStereogenic),
            Some(true) | None => assign_double_bond_descriptor(mol, id, stereo, options),
        },
        StereoElementKind::Axis(stereo) => assign_axis_descriptor(mol, id, stereo, options),
    };
    match assignment {
        Ok(descriptor) => CipElementAssignment::Assigned(descriptor),
        Err(CipAssignmentIssue::UnresolvedPriority { .. }) => CipElementAssignment::Deferred,
        Err(issue) => CipElementAssignment::Issue(issue),
    }
}

fn double_bond_cip_stereogenic(mol: &Molecule, stereo: &DoubleBondStereo) -> Option<bool> {
    let bond = mol.bond(stereo.bond).ok()?;
    if bond.order != BondOrder::Double {
        return None;
    }
    if mol.bond_is_aromatic(stereo.bond).ok().flatten() == Some(true)
        || double_bond_between_aromatic_atoms(mol, bond)
    {
        return Some(false);
    }
    if bond_in_ring_smaller_than(mol, stereo.bond, 8) {
        return Some(false);
    }
    if double_bond_is_in_ring(mol, stereo.bond) && double_bond_has_noncarbon_endpoint(mol, bond) {
        return Some(false);
    }
    Some(true)
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
    assign_tetrahedral_descriptor_with_deferred_rule6(mol, element, stereo, options, false)
}

fn assign_tetrahedral_descriptor_with_deferred_rule6(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &TetrahedralStereo,
    options: CipAssignmentOptions,
    allow_single_ring_tied_pair_rule6: bool,
) -> CipResult<StereoDescriptor> {
    let orientation = stereo
        .orientation
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let ranked = ranked_tetrahedral_carriers(
        mol,
        element,
        stereo.center,
        &stereo.carriers,
        orientation,
        options,
        allow_single_ring_tied_pair_rule6,
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
    if bond_in_ring_smaller_than(mol, stereo.bond, 8) {
        return Err(CipAssignmentIssue::UnresolvedPriority { element });
    }
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
    let left_top = left_ranked
        .carriers
        .first()
        .copied()
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let right_ranked = ranked_carriers(
        mol,
        element,
        stereo.right,
        &right_carriers,
        options,
        true,
        false,
    )?;
    let right_top = right_ranked
        .carriers
        .first()
        .copied()
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;

    let mut top_relation = orientation;
    if stereo.left_carrier != left_top {
        top_relation = top_relation.inverted();
    }
    if stereo.right_carrier != right_top {
        top_relation = top_relation.inverted();
    }
    let pseudo_sequence =
        left_ranked.pseudo_asymmetric_ordering != right_ranked.pseudo_asymmetric_ordering;
    Ok(match (top_relation, pseudo_sequence) {
        (DoubleBondOrientation::Together, true) => StereoDescriptor::SeqCis,
        (DoubleBondOrientation::Opposite, true) => StereoDescriptor::SeqTrans,
        (DoubleBondOrientation::Together, false) => StereoDescriptor::Z,
        (DoubleBondOrientation::Opposite, false) => StereoDescriptor::E,
    })
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
    let left_top = left_ranked
        .carriers
        .first()
        .copied()
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let right_ranked = ranked_carriers(
        mol,
        element,
        right,
        &axis_endpoint_carriers(mol, right, left, stereo.axis),
        options,
        true,
        true,
    )?;
    let right_top = right_ranked
        .carriers
        .first()
        .copied()
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let mut top_orientation = orientation;
    if left_reference != left_top {
        top_orientation = top_orientation.inverted();
    }
    if right_reference != right_top {
        top_orientation = top_orientation.inverted();
    }
    let pseudo_axis =
        left_ranked.pseudo_asymmetric_ordering || right_ranked.pseudo_asymmetric_ordering;
    Ok(match (top_orientation, pseudo_axis) {
        (AxisOrientation::CounterClockwise, true) => StereoDescriptor::LowerM,
        (AxisOrientation::Clockwise, true) => StereoDescriptor::LowerP,
        (AxisOrientation::CounterClockwise, false) => StereoDescriptor::M,
        (AxisOrientation::Clockwise, false) => StereoDescriptor::P,
    })
}

fn axis_reference_carriers(
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

fn axis_endpoint_carriers(
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
    normalize_all_carbon_aromatic: bool,
) -> CipResult<RankedCarriers> {
    let signatures = carrier_signatures(
        mol,
        element,
        root,
        carriers,
        options,
        allow_auxiliary_descriptors,
        normalize_all_carbon_aromatic,
    )?;
    rank_carrier_signatures(element, &signatures, None)
}

fn ranked_tetrahedral_carriers(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    carriers: &[StereoCarrier],
    orientation: TetrahedralOrientation,
    options: CipAssignmentOptions,
    allow_single_ring_tied_pair_rule6: bool,
) -> CipResult<RankedCarriers> {
    let signatures = carrier_signatures(
        mol,
        element,
        root,
        carriers,
        options,
        allow_single_ring_tied_pair_rule6,
        false,
    )?;
    match rank_carrier_signatures(element, &signatures, None) {
        Ok(ranked) => Ok(ranked),
        Err(CipAssignmentIssue::UnresolvedPriority { .. }) if carriers.len() == 4 => {
            rank_tetrahedral_signatures_with_rule6(
                mol,
                element,
                root,
                &signatures,
                orientation,
                allow_single_ring_tied_pair_rule6,
            )
        }
        Err(issue) => Err(issue),
    }
}

pub(super) fn assign_deferred_tetrahedral_rule6(
    mol: &Molecule,
    pending: &[(StereoElementId, StereoElement)],
    options: CipAssignmentOptions,
) -> CipResult<Vec<(StereoElementId, StereoDescriptor)>> {
    let mut assignments = Vec::new();
    for (id, element) in pending {
        let StereoElementKind::Tetrahedral(stereo) = &element.kind else {
            continue;
        };
        if !element.is_specified() {
            continue;
        }
        match assign_tetrahedral_descriptor_with_deferred_rule6(mol, *id, stereo, options, true) {
            Ok(descriptor) => assignments.push((*id, descriptor)),
            Err(CipAssignmentIssue::UnresolvedPriority { .. }) => {}
            Err(issue) => return Err(issue),
        }
    }
    Ok(assignments)
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
    let orientation = stereo
        .orientation
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let options = complete_final_tie_options(mol, stereo.center, options);
    let signatures = carrier_signatures(
        mol,
        element,
        stereo.center,
        &stereo.carriers,
        options,
        true,
        false,
    )?;
    match rank_tetrahedral_signatures_with_rule6(
        mol,
        element,
        stereo.center,
        &signatures,
        orientation,
        true,
    ) {
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
    normalize_all_carbon_aromatic: bool,
) -> CipResult<bool> {
    if carriers.len() < 2 {
        return Ok(false);
    }
    let options = complete_final_tie_options(mol, root, options);
    let signatures = carrier_signatures(
        mol,
        element,
        root,
        carriers,
        options,
        true,
        normalize_all_carbon_aromatic,
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

fn complete_final_tie_options(
    mol: &Molecule,
    root: AtomId,
    mut options: CipAssignmentOptions,
) -> CipAssignmentOptions {
    options.max_depth = options.max_depth.max(connected_atom_count(mol, root));
    options
}

fn connected_atom_count(mol: &Molecule, root: AtomId) -> usize {
    if mol.atom(root).is_err() {
        return 0;
    }
    let mut seen = HashSet::new();
    let mut queue = VecDeque::from([root]);
    while let Some(atom) = queue.pop_front() {
        if !seen.insert(atom) {
            continue;
        }
        if let Ok(incident) = mol.incident_bonds(atom) {
            for (_, bond) in incident {
                let neighbor = bond.other_atom(atom);
                if !seen.contains(&neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }
    }
    seen.len()
}

fn carrier_signatures(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    carriers: &[StereoCarrier],
    options: CipAssignmentOptions,
    allow_auxiliary_descriptors: bool,
    normalize_all_carbon_aromatic: bool,
) -> CipResult<Vec<(StereoCarrier, LigandSignature)>> {
    let cip_bond_orders = CipBondOrders::new(mol, normalize_all_carbon_aromatic);
    let atomic_number_fractions = cip_atomic_number_fractions(mol, &cip_bond_orders);
    if allow_auxiliary_descriptors {
        let descriptor_context = DescriptorContext::new(element, AuxiliaryDescriptorMode::Collect);
        let aux_graph = build_auxiliary_graph(
            mol,
            element,
            root,
            options,
            &atomic_number_fractions,
            &cip_bond_orders,
        )?;
        collect_auxiliary_occurrences_from_molecule(mol, &descriptor_context, &aux_graph);
        precompute_auxiliary_descriptors(
            mol,
            &descriptor_context,
            &aux_graph,
            options,
            &atomic_number_fractions,
            &cip_bond_orders,
        );
        let descriptor_context = descriptor_context.with_mode(AuxiliaryDescriptorMode::Precomputed);
        let build_context = LigandBuildContext {
            mol,
            element,
            descriptor_context: &descriptor_context,
            options,
            atomic_number_fractions: &atomic_number_fractions,
            cip_bond_orders: &cip_bond_orders,
        };
        let signatures = build_carrier_signatures(&build_context, root, carriers)?;
        return Ok(signatures);
    }
    let descriptor_context = DescriptorContext::new(element, AuxiliaryDescriptorMode::Disabled);
    let build_context = LigandBuildContext {
        mol,
        element,
        descriptor_context: &descriptor_context,
        options,
        atomic_number_fractions: &atomic_number_fractions,
        cip_bond_orders: &cip_bond_orders,
    };
    build_carrier_signatures(&build_context, root, carriers)
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
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    signatures: &[(StereoCarrier, LigandSignature)],
    orientation: TetrahedralOrientation,
    allow_single_ring_tied_pair_rule6: bool,
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
        _ if allow_single_ring_tied_pair_rule6 => rank_single_ring_tied_pair_with_rule6(
            mol,
            element,
            root,
            orientation,
            signatures,
            &groups,
        ),
        _ => Err(CipAssignmentIssue::UnresolvedPriority { element }),
    }
}

fn rank_single_ring_tied_pair_with_rule6(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    orientation: TetrahedralOrientation,
    signatures: &[(StereoCarrier, LigandSignature)],
    groups: &[Vec<usize>],
) -> CipResult<RankedCarriers> {
    let tied_groups = groups
        .iter()
        .filter(|group| group.len() > 1)
        .collect::<Vec<_>>();
    if tied_groups.len() != 1 || tied_groups[0].len() != 2 {
        return Err(CipAssignmentIssue::UnresolvedPriority { element });
    }
    let left = carrier_rule6_atom(signatures[tied_groups[0][0]].0)
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let right = carrier_rule6_atom(signatures[tied_groups[0][1]].0)
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let Some(path) = shortest_path_excluding_root(mol, left, right, root) else {
        return Err(CipAssignmentIssue::UnresolvedPriority { element });
    };
    let path_length = path.len().saturating_sub(1);
    let tied_pair_descriptor_class = tied_groups[0]
        .iter()
        .filter_map(|index| tree_descriptor_class(&signatures[*index].1.root))
        .max();
    let outside_tied_pair_descriptor_class = groups
        .iter()
        .filter(|group| !std::ptr::eq(*group, tied_groups[0]))
        .flatten()
        .filter_map(|index| tree_descriptor_class(&signatures[*index].1.root))
        .max();
    let tied_pair_descriptor_refs_match =
        descriptor_ref_counts(&signatures[tied_groups[0][0]].1.root)
            == descriptor_ref_counts(&signatures[tied_groups[0][1]].1.root);
    let reference = if tied_pair_descriptor_class.is_some() {
        if tied_pair_descriptor_refs_match {
            return Err(CipAssignmentIssue::UnresolvedPriority { element });
        }
        if left.raw() >= right.raw() {
            left
        } else {
            right
        }
    } else {
        match outside_tied_pair_descriptor_class {
            Some(DescriptorClass::Absolute) => {
                if left.raw() >= right.raw() {
                    left
                } else {
                    right
                }
            }
            Some(DescriptorClass::Pseudo) => {
                if left.raw() <= right.raw() {
                    left
                } else {
                    right
                }
            }
            None if path_length == 2 => {
                match path
                    .get(1)
                    .and_then(|center| tetrahedral_orientation_for_center(mol, *center))
                {
                    Some(other_orientation) if other_orientation != orientation => right,
                    _ => left,
                }
            }
            None if mol.stereo_elements().all(|(id, _)| id == element)
                && ring_path_is_unsubstituted_bridge(mol, &path, root) =>
            {
                return Err(CipAssignmentIssue::UnresolvedPriority { element });
            }
            None if left.raw() >= right.raw() => left,
            None => right,
        }
    };
    let mut ranked = rank_carrier_signatures(element, signatures, Some(reference))?;
    ranked.pseudo_asymmetric_ordering = !matches!(
        (
            tied_pair_descriptor_class,
            outside_tied_pair_descriptor_class,
        ),
        (Some(DescriptorClass::Absolute), _) | (None, Some(DescriptorClass::Absolute))
    );
    Ok(ranked)
}

fn descriptor_ref_counts(tree: &LigandTree) -> (usize, usize) {
    let own = match tree.priority.descriptor.and_then(descriptor_ref) {
        Some(DescriptorRef::R) => (1, 0),
        Some(DescriptorRef::S) => (0, 1),
        None => (0, 0),
    };
    tree.children
        .iter()
        .map(descriptor_ref_counts)
        .fold(own, |left, right| (left.0 + right.0, left.1 + right.1))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum DescriptorClass {
    Pseudo,
    Absolute,
}

fn tree_descriptor_class(tree: &LigandTree) -> Option<DescriptorClass> {
    let own = tree.priority.descriptor.and_then(descriptor_class);
    tree.children
        .iter()
        .filter_map(tree_descriptor_class)
        .chain(own)
        .max()
}

fn descriptor_class(descriptor: StereoDescriptor) -> Option<DescriptorClass> {
    match descriptor {
        StereoDescriptor::R
        | StereoDescriptor::S
        | StereoDescriptor::M
        | StereoDescriptor::P
        | StereoDescriptor::SeqCis
        | StereoDescriptor::SeqTrans => Some(DescriptorClass::Absolute),
        StereoDescriptor::LowerR
        | StereoDescriptor::LowerS
        | StereoDescriptor::LowerM
        | StereoDescriptor::LowerP => Some(DescriptorClass::Pseudo),
        StereoDescriptor::E | StereoDescriptor::Z => None,
    }
}

fn tetrahedral_orientation_for_center(
    mol: &Molecule,
    center: AtomId,
) -> Option<TetrahedralOrientation> {
    mol.stereo_elements()
        .find_map(|(_, element)| match &element.kind {
            StereoElementKind::Tetrahedral(stereo) if stereo.center == center => stereo.orientation,
            _ => None,
        })
}

fn shortest_path_excluding_root(
    mol: &Molecule,
    left: AtomId,
    right: AtomId,
    root: AtomId,
) -> Option<Vec<AtomId>> {
    let mut seen = Vec::new();
    let mut queue = VecDeque::from([(left, vec![left])]);
    while let Some((atom, path)) = queue.pop_front() {
        if atom == right {
            return Some(path);
        }
        if atom == root || seen.contains(&atom) {
            continue;
        }
        seen.push(atom);
        if let Ok(incident) = mol.incident_bonds(atom) {
            for (_, bond) in incident {
                let neighbor = bond.other_atom(atom);
                if neighbor != root && !seen.contains(&neighbor) {
                    let mut next_path = path.clone();
                    next_path.push(neighbor);
                    queue.push_back((neighbor, next_path));
                }
            }
        }
    }
    None
}

fn ring_path_is_unsubstituted_bridge(mol: &Molecule, path: &[AtomId], root: AtomId) -> bool {
    path.iter().all(|atom| {
        mol.incident_bonds(*atom)
            .map(|incident| {
                incident.into_iter().all(|(_, bond)| {
                    let neighbor = bond.other_atom(*atom);
                    neighbor == root || path.contains(&neighbor)
                })
            })
            .unwrap_or(false)
    })
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
