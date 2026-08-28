use crate::*;

use super::chemistry::*;
use super::descriptors::{rdkit_default_atom_index, rdkit_default_bond_count};
use super::io::{IndexedSmallRecord, IndexedStereoPerceptionRecord};
use super::smiles::offset_object_u64;

pub(crate) fn stereo_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = &record.molecule;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "stereo_elements": stereo_elements_json(mol),
        "stereo_groups": stereo_groups_json(mol),
    })
}

pub(crate) fn stereo_perception_record_json(record: &mut IndexedSmallRecord) -> Value {
    let source_stereo_elements = record
        .molecule
        .stereo_elements()
        .map(|(id, _)| id)
        .collect::<Vec<_>>();
    if record.molecule.perceive().is_err() {
        let mol = &record.molecule;
        return json!({
            "record_index": record.record_index,
            "status": "perception_error",
            "title": record.title,
            "atom_count": mol.atom_count(),
            "bond_count": mol.bond_count(),
        });
    }
    let candidates = stereo::detect_stereo_candidates(&record.molecule);
    let positions = kekule::structure::Positions::zeros(record.molecule.atom_count());
    let mut editor = record.molecule.edit();
    let result = stereo::materialize_coordinate_stereo(&mut editor, &positions);
    if result.is_ok() {
        record.molecule = editor
            .finish()
            .expect("coordinate stereo materialization preserves publication invariants");
    }
    let mol = &record.molecule;
    match result {
        Ok(report) => json!({
            "record_index": record.record_index,
            "status": "ok",
            "title": record.title,
            "atom_count": mol.atom_count(),
            "bond_count": mol.bond_count(),
            "report": stereo_perception_benchmark_report_json(
                mol,
                &source_stereo_elements,
                &candidates,
                &report,
            ),
            "stereo_elements": stereo_elements_json(mol),
            "stereo_groups": stereo_groups_json(mol),
        }),
        Err(error) => json!({
            "record_index": record.record_index,
            "status": "perception_error",
            "title": record.title,
            "atom_count": mol.atom_count(),
            "bond_count": mol.bond_count(),
            "candidates": candidates.iter().map(stereo_candidate_json).collect::<Vec<_>>(),
            "source_stereo_element_indices": source_stereo_elements
                .iter()
                .map(|id| id.raw())
                .collect::<Vec<_>>(),
            "error": coordinate_stereo_error_json(&error),
            "stereo_elements": stereo_elements_json(mol),
            "stereo_groups": stereo_groups_json(mol),
        }),
    }
}

pub(crate) fn stereo_perception_group_record_json(
    record: &mut IndexedStereoPerceptionRecord,
) -> Value {
    if record.components.is_empty() {
        return json!({
            "record_index": record.record_index,
            "status": "parse_error",
            "title": record.title,
            "atom_count": 0,
            "bond_count": 0,
        });
    }

    let mut atom_count = 0u64;
    let mut bond_count = 0u64;
    let mut element_count = 0u64;
    let mut group_count = 0u64;
    let mut assembled_count = 0u64;
    let mut assembled_elements = Vec::new();
    let mut candidates = Vec::new();
    let mut created_element_indices = Vec::new();
    let mut issues = Vec::new();
    let mut stereo_elements = Vec::new();
    let mut stereo_groups = Vec::new();

    for component in &record.components {
        let mut component_record = IndexedSmallRecord {
            record_index: record.record_index,
            title: record.title.clone(),
            molecule: component.clone(),
            sdf_fields: BTreeMap::new(),
        };
        let mut value = stereo_perception_record_json(&mut component_record);
        if value.get("status").and_then(Value::as_str) != Some("ok") {
            return json!({
                "record_index": record.record_index,
                "status": value.get("status").cloned().unwrap_or_else(|| json!("perception_error")),
                "title": record.title,
                "atom_count": record.components.iter().map(|molecule| molecule.atom_count()).sum::<usize>(),
                "bond_count": record.components.iter().map(|molecule| molecule.bond_count()).sum::<usize>(),
            });
        }
        let object = value
            .as_object_mut()
            .expect("stereo record must be an object");
        let component_atom_count = object
            .get("atom_count")
            .and_then(Value::as_u64)
            .unwrap_or_default();
        let component_bond_count = object
            .get("bond_count")
            .and_then(Value::as_u64)
            .unwrap_or_default();
        let mut report = match object.remove("report") {
            Some(Value::Object(report)) => report,
            _ => panic!("successful stereo record must contain a report"),
        };

        let mut component_assembled = take_array(&mut report, "assembled_elements");
        for element in &mut component_assembled {
            offset_object_u64(element, "index", assembled_count);
            offset_stereo_references(element, atom_count, bond_count, element_count, group_count);
        }
        assembled_count += component_assembled.len() as u64;
        assembled_elements.extend(component_assembled);

        let mut component_candidates = take_array(&mut report, "candidates");
        for candidate in &mut component_candidates {
            offset_stereo_references(
                candidate,
                atom_count,
                bond_count,
                element_count,
                group_count,
            );
        }
        candidates.extend(component_candidates);

        for created in take_array(&mut report, "created_element_indices") {
            if let Some(index) = created.as_u64() {
                created_element_indices.push(json!(index + element_count));
            }
        }
        let mut component_issues = take_array(&mut report, "issues");
        for issue in &mut component_issues {
            offset_stereo_references(issue, atom_count, bond_count, element_count, group_count);
        }
        issues.extend(component_issues);

        let mut component_elements = take_array(object, "stereo_elements");
        for element in &mut component_elements {
            offset_object_u64(element, "index", element_count);
            offset_stereo_references(element, atom_count, bond_count, element_count, group_count);
        }
        let component_element_count = component_elements.len() as u64;
        stereo_elements.extend(component_elements);

        let mut component_groups = take_array(object, "stereo_groups");
        for group in &mut component_groups {
            offset_object_u64(group, "index", group_count);
            if let Some(members) = group
                .as_object_mut()
                .and_then(|object| object.get_mut("members"))
                .and_then(Value::as_array_mut)
            {
                for member in members {
                    if let Some(index) = member.as_u64() {
                        *member = json!(index + element_count);
                    }
                }
            }
        }
        let component_group_count = component_groups.len() as u64;
        stereo_groups.extend(component_groups);

        atom_count += component_atom_count;
        bond_count += component_bond_count;
        element_count += component_element_count;
        group_count += component_group_count;
    }

    candidates.sort_by_key(
        |candidate| match candidate.get("type").and_then(Value::as_str) {
            Some("tetrahedral") => (
                0u8,
                candidate
                    .get("center_atom_index")
                    .and_then(Value::as_u64)
                    .unwrap_or(u64::MAX),
            ),
            Some("double_bond") => (
                1u8,
                candidate
                    .get("center_bond_index")
                    .and_then(Value::as_u64)
                    .unwrap_or(u64::MAX),
            ),
            _ => (u8::MAX, u64::MAX),
        },
    );

    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": atom_count,
        "bond_count": bond_count,
        "report": {
            "is_ok": issues.is_empty(),
            "candidates": candidates,
            "issues": issues,
            "assembled_elements": assembled_elements,
            "created_element_indices": created_element_indices,
        },
        "stereo_elements": stereo_elements,
        "stereo_groups": stereo_groups,
    })
}

pub(crate) fn stereo_cip_record_json(
    record: &mut IndexedSmallRecord,
    remove_plain_hydrogens: bool,
) -> Option<Value> {
    if record.molecule.perceive().is_err() {
        return None;
    }
    if stereo::validate_stereo(&record.molecule).is_err() {
        return None;
    }
    stereo::assign_cip_descriptors(&mut record.molecule).ok()?;
    let mol = &record.molecule;
    let atom_index = rdkit_default_atom_index(mol, remove_plain_hydrogens);
    let atom_descriptors = cip_atom_descriptors_json(mol, &atom_index);
    let bond_descriptors = cip_bond_descriptors_json(mol, &atom_index);
    if atom_descriptors.is_empty() && bond_descriptors.is_empty() {
        return None;
    }
    Some(json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": atom_index.len(),
        "bond_count": rdkit_default_bond_count(mol, &atom_index),
        "atom_descriptors": atom_descriptors,
        "bond_descriptors": bond_descriptors,
    }))
}

pub(crate) fn cip_atom_descriptors_json(
    mol: &Molecule,
    atom_index: &BTreeMap<AtomId, u64>,
) -> Vec<Value> {
    let mut descriptors = mol
        .stereo_elements()
        .filter_map(|(id, element)| match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => mol
                .cip_descriptor(id)
                .ok()
                .flatten()
                .and_then(|descriptor| {
                    let atom_index = *atom_index.get(&stereo.center)?;
                    Some(json!({
                        "atom_index": atom_index,
                        "descriptor": stereo_descriptor_json(descriptor),
                    }))
                }),
            StereoElementKind::Axis(_) | StereoElementKind::DoubleBond(_) => None,
        })
        .collect::<Vec<_>>();
    descriptors.sort_by_key(|value| {
        value
            .get("atom_index")
            .and_then(Value::as_u64)
            .unwrap_or(u64::MAX)
    });
    descriptors
}

pub(crate) fn cip_bond_descriptors_json(
    mol: &Molecule,
    atom_index: &BTreeMap<AtomId, u64>,
) -> Vec<Value> {
    let mut descriptors = mol
        .stereo_elements()
        .filter_map(|(id, element)| match &element.kind {
            StereoElementKind::DoubleBond(stereo) => mol
                .cip_descriptor(id)
                .ok()
                .flatten()
                .and_then(|descriptor| {
                    let begin_atom_index = *atom_index.get(&stereo.left)?;
                    let end_atom_index = *atom_index.get(&stereo.right)?;
                    Some(json!({
                        "begin_atom_index": begin_atom_index,
                        "end_atom_index": end_atom_index,
                        "descriptor": stereo_descriptor_json(descriptor),
                    }))
                }),
            StereoElementKind::Axis(stereo) => {
                mol.cip_descriptor(id)
                    .ok()
                    .flatten()
                    .and_then(|descriptor| {
                        let bond = mol.bond(stereo.axis).ok()?;
                        let (begin, end) = bond.endpoints();
                        let begin_atom_index = *atom_index.get(&begin)?;
                        let end_atom_index = *atom_index.get(&end)?;
                        Some(json!({
                            "begin_atom_index": begin_atom_index,
                            "end_atom_index": end_atom_index,
                            "descriptor": stereo_descriptor_json(descriptor),
                        }))
                    })
            }
            StereoElementKind::Tetrahedral(_) => None,
        })
        .collect::<Vec<_>>();
    descriptors.sort_by(|left, right| {
        let left_key = (
            left.get("begin_atom_index")
                .and_then(Value::as_u64)
                .unwrap_or(u64::MAX),
            left.get("end_atom_index")
                .and_then(Value::as_u64)
                .unwrap_or(u64::MAX),
        );
        let right_key = (
            right
                .get("begin_atom_index")
                .and_then(Value::as_u64)
                .unwrap_or(u64::MAX),
            right
                .get("end_atom_index")
                .and_then(Value::as_u64)
                .unwrap_or(u64::MAX),
        );
        left_key.cmp(&right_key).then_with(|| {
            left.get("descriptor")
                .and_then(Value::as_str)
                .unwrap_or("")
                .cmp(
                    right
                        .get("descriptor")
                        .and_then(Value::as_str)
                        .unwrap_or(""),
                )
        })
    });
    descriptors
}

pub(crate) fn stereo_perception_benchmark_report_json(
    mol: &Molecule,
    source_stereo_elements: &[StereoElementId],
    candidates: &[StereoCandidate],
    report: &CoordinateStereoMaterializationReport,
) -> Value {
    let assembled_elements = source_stereo_elements
        .iter()
        .filter_map(|id| {
            mol.stereo_element(*id)
                .ok()
                .map(|element| stereo_element_json(u64::from(id.raw()), element, None))
        })
        .collect::<Vec<_>>();
    json!({
        "is_ok": true,
        "candidates": candidates.iter().map(stereo_candidate_json).collect::<Vec<_>>(),
        "issues": [],
        "assembled_elements": assembled_elements,
        "created_element_indices": report
            .created_elements
            .iter()
            .map(|id| id.raw())
            .collect::<Vec<_>>(),
    })
}

fn take_array(object: &mut serde_json::Map<String, Value>, key: &str) -> Vec<Value> {
    match object.remove(key) {
        Some(Value::Array(values)) => values,
        _ => Vec::new(),
    }
}

fn offset_stereo_references(
    value: &mut Value,
    atom_offset: u64,
    bond_offset: u64,
    element_offset: u64,
    group_offset: u64,
) {
    match value {
        Value::Array(values) => {
            for value in values {
                offset_stereo_references(
                    value,
                    atom_offset,
                    bond_offset,
                    element_offset,
                    group_offset,
                );
            }
        }
        Value::Object(object) => {
            for (key, value) in object {
                let offset = if key == "atom_index" || key.ends_with("_atom_index") {
                    atom_offset
                } else if key == "bond_index" || key.ends_with("_bond_index") {
                    bond_offset
                } else if key == "element_index" || key.ends_with("_element_index") {
                    element_offset
                } else if key == "group_index" || key.ends_with("_group_index") {
                    group_offset
                } else {
                    offset_stereo_references(
                        value,
                        atom_offset,
                        bond_offset,
                        element_offset,
                        group_offset,
                    );
                    continue;
                };
                if let Some(index) = value.as_u64() {
                    *value = json!(index + offset);
                }
            }
        }
        _ => {}
    }
}

pub(crate) fn stereo_candidate_json(candidate: &StereoCandidate) -> Value {
    match candidate {
        StereoCandidate::Tetrahedral { center, carriers } => json!({
            "type": "tetrahedral",
            "center_atom_index": center.raw(),
            "carriers": carriers.iter().map(stereo_carrier_json).collect::<Vec<_>>(),
        }),
        StereoCandidate::DoubleBond {
            bond,
            left,
            right,
            left_carriers,
            right_carriers,
        } => json!({
            "type": "double_bond",
            "center_bond_index": bond.raw(),
            "left_atom_index": left.raw(),
            "right_atom_index": right.raw(),
            "left_carriers": left_carriers.iter().map(stereo_carrier_json).collect::<Vec<_>>(),
            "right_carriers": right_carriers.iter().map(stereo_carrier_json).collect::<Vec<_>>(),
        }),
    }
}

pub(crate) fn coordinate_stereo_error_json(error: &CoordinateStereoError) -> Value {
    let issues = match error {
        CoordinateStereoError::InvalidStereo(error) => error
            .issues
            .iter()
            .map(|issue| {
                json!({
                    "type": "invalid_stereo",
                    "issue": stereo_validation_issue_json(issue),
                })
            })
            .collect::<Vec<_>>(),
        CoordinateStereoError::CouldNotCreateElement(error) => vec![json!({
            "type": "could_not_create_element",
            "error": format!("{error:?}"),
        })],
        _ => vec![json!({
            "type": "coordinate_stereo_error",
            "error": format!("{error:?}"),
        })],
    };
    json!({ "issues": issues })
}

pub(crate) fn stereo_validation_issue_json(issue: &StereoValidationIssue) -> Value {
    match issue {
        StereoValidationIssue::MissingStereoAtom { element, atom } => json!({
            "type": "missing_stereo_atom",
            "element_index": element.raw(),
            "atom_index": atom.raw(),
        }),
        StereoValidationIssue::MissingStereoBond { element, bond } => json!({
            "type": "missing_stereo_bond",
            "element_index": element.raw(),
            "bond_index": bond.raw(),
        }),
        StereoValidationIssue::InvalidTetrahedralCarrierCount {
            element,
            center,
            carrier_count,
        } => json!({
            "type": "invalid_tetrahedral_carrier_count",
            "element_index": element.raw(),
            "center_atom_index": center.raw(),
            "carrier_count": carrier_count,
        }),
        StereoValidationIssue::DuplicateTetrahedralCarrier {
            element,
            center,
            carrier,
        } => json!({
            "type": "duplicate_tetrahedral_carrier",
            "element_index": element.raw(),
            "center_atom_index": center.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
        StereoValidationIssue::TetrahedralCarrierNotAdjacent {
            element,
            center,
            carrier,
        } => json!({
            "type": "tetrahedral_carrier_not_adjacent",
            "element_index": element.raw(),
            "center_atom_index": center.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
        StereoValidationIssue::InvalidDoubleBondOrder {
            element,
            bond,
            order,
        } => json!({
            "type": "invalid_double_bond_order",
            "element_index": element.raw(),
            "bond_index": bond.raw(),
            "bond_order": bond_order_json(*order),
        }),
        StereoValidationIssue::DoubleBondFocusMismatch {
            element,
            bond,
            left,
            right,
        } => json!({
            "type": "double_bond_focus_mismatch",
            "element_index": element.raw(),
            "bond_index": bond.raw(),
            "left_atom_index": left.raw(),
            "right_atom_index": right.raw(),
        }),
        StereoValidationIssue::DoubleBondCarrierIsFocusAtom {
            element,
            endpoint,
            carrier,
        } => json!({
            "type": "double_bond_carrier_is_focus_atom",
            "element_index": element.raw(),
            "endpoint_atom_index": endpoint.raw(),
            "carrier_atom_index": carrier.raw(),
        }),
        StereoValidationIssue::DoubleBondCarrierNotAdjacent {
            element,
            endpoint,
            carrier,
        } => json!({
            "type": "double_bond_carrier_not_adjacent",
            "element_index": element.raw(),
            "endpoint_atom_index": endpoint.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
        StereoValidationIssue::UnsupportedDoubleBondCarrier {
            element,
            endpoint,
            carrier,
        } => json!({
            "type": "unsupported_double_bond_carrier",
            "element_index": element.raw(),
            "endpoint_atom_index": endpoint.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
        StereoValidationIssue::InvalidAxisCarrierCount {
            element,
            axis,
            carrier_count,
        } => json!({
            "type": "invalid_axis_carrier_count",
            "element_index": element.raw(),
            "axis_bond_index": axis.raw(),
            "carrier_count": carrier_count,
        }),
        StereoValidationIssue::AxisCarrierIsFocusAtom {
            element,
            axis,
            carrier,
        } => json!({
            "type": "axis_carrier_is_focus_atom",
            "element_index": element.raw(),
            "axis_bond_index": axis.raw(),
            "carrier_atom_index": carrier.raw(),
        }),
        StereoValidationIssue::AxisCarrierNotAdjacent {
            element,
            axis,
            carrier,
        } => json!({
            "type": "axis_carrier_not_adjacent",
            "element_index": element.raw(),
            "axis_bond_index": axis.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
        StereoValidationIssue::UnsupportedAxisCarrier {
            element,
            axis,
            carrier,
        } => json!({
            "type": "unsupported_axis_carrier",
            "element_index": element.raw(),
            "axis_bond_index": axis.raw(),
            "carrier": stereo_carrier_json(carrier),
        }),
    }
}

pub(crate) fn stereo_elements_json(mol: &Molecule) -> Vec<Value> {
    mol.stereo_elements()
        .map(|(id, element)| {
            stereo_element_json(
                u64::from(id.raw()),
                element,
                mol.cip_descriptor(id).ok().flatten(),
            )
        })
        .collect()
}

pub(crate) fn stereo_element_json(
    index: u64,
    element: &StereoElement,
    descriptor: Option<StereoDescriptor>,
) -> Value {
    let mut object = serde_json::Map::new();
    object.insert("index".to_owned(), json!(index));
    if let Some(group) = element.group {
        object.insert("group_index".to_owned(), json!(group.raw()));
    }
    if let Some(descriptor) = descriptor {
        object.insert(
            "descriptor".to_owned(),
            json!(stereo_descriptor_json(descriptor)),
        );
    }
    match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            object.insert("type".to_owned(), json!("tetrahedral"));
            object.insert("center_atom_index".to_owned(), json!(stereo.center.raw()));
            object.insert(
                "carriers".to_owned(),
                Value::Array(
                    stereo
                        .carriers
                        .iter()
                        .map(stereo_carrier_json)
                        .collect::<Vec<_>>(),
                ),
            );
            object.insert(
                "orientation".to_owned(),
                json!(stereo.orientation.map(tetrahedral_orientation_json)),
            );
        }
        StereoElementKind::DoubleBond(stereo) => {
            object.insert("type".to_owned(), json!("double_bond"));
            object.insert("center_bond_index".to_owned(), json!(stereo.bond.raw()));
            object.insert("left_atom_index".to_owned(), json!(stereo.left.raw()));
            object.insert("right_atom_index".to_owned(), json!(stereo.right.raw()));
            object.insert(
                "left_carrier".to_owned(),
                stereo_carrier_json(&stereo.left_carrier),
            );
            object.insert(
                "right_carrier".to_owned(),
                stereo_carrier_json(&stereo.right_carrier),
            );
            object.insert(
                "orientation".to_owned(),
                json!(stereo.orientation.map(double_bond_orientation_json)),
            );
        }
        StereoElementKind::Axis(stereo) => {
            object.insert("type".to_owned(), json!("axis"));
            object.insert("axis_bond_index".to_owned(), json!(stereo.axis.raw()));
            object.insert(
                "carriers".to_owned(),
                Value::Array(
                    stereo
                        .carriers
                        .iter()
                        .map(stereo_carrier_json)
                        .collect::<Vec<_>>(),
                ),
            );
            object.insert(
                "orientation".to_owned(),
                json!(stereo.orientation.map(axis_orientation_json)),
            );
        }
    }
    Value::Object(object)
}

pub(crate) fn stereo_groups_json(mol: &Molecule) -> Vec<Value> {
    mol.stereo_groups()
        .map(|(id, group)| stereo_group_json(id.raw(), group))
        .collect()
}

pub(crate) fn stereo_group_json(index: u32, group: &StereoGroup) -> Value {
    json!({
        "index": index,
        "kind": stereo_group_kind_json(group.kind),
        "members": group.members.iter().map(|member| member.raw()).collect::<Vec<_>>(),
    })
}

pub(crate) fn stereo_carrier_json(carrier: &StereoCarrier) -> Value {
    match carrier {
        StereoCarrier::Atom(atom) => json!({ "atom_index": atom.raw() }),
        StereoCarrier::ImplicitHydrogen => json!({ "implicit_hydrogen": true }),
        StereoCarrier::ImplicitLonePair => json!({ "implicit_lone_pair": true }),
    }
}

pub(crate) fn stereo_descriptor_json(descriptor: StereoDescriptor) -> &'static str {
    match descriptor {
        StereoDescriptor::R => "R",
        StereoDescriptor::S => "S",
        StereoDescriptor::LowerR => "r",
        StereoDescriptor::LowerS => "s",
        StereoDescriptor::SeqTrans => "seqTrans",
        StereoDescriptor::SeqCis => "seqCis",
        StereoDescriptor::E => "E",
        StereoDescriptor::Z => "Z",
        StereoDescriptor::M => "M",
        StereoDescriptor::P => "P",
        StereoDescriptor::LowerM => "m",
        StereoDescriptor::LowerP => "p",
    }
}

pub(crate) fn stereo_group_kind_json(kind: StereoGroupKind) -> &'static str {
    match kind {
        StereoGroupKind::Absolute => "absolute",
        StereoGroupKind::Relative => "relative",
        StereoGroupKind::Racemic => "racemic",
        StereoGroupKind::And => "and",
        StereoGroupKind::Or => "or",
    }
}

pub(crate) fn tetrahedral_orientation_json(orientation: TetrahedralOrientation) -> &'static str {
    match orientation {
        TetrahedralOrientation::Clockwise => "clockwise",
        TetrahedralOrientation::CounterClockwise => "counter_clockwise",
    }
}

pub(crate) fn double_bond_orientation_json(orientation: DoubleBondOrientation) -> &'static str {
    match orientation {
        DoubleBondOrientation::Together => "together",
        DoubleBondOrientation::Opposite => "opposite",
    }
}

pub(crate) fn axis_orientation_json(orientation: AxisOrientation) -> &'static str {
    match orientation {
        AxisOrientation::Clockwise => "clockwise",
        AxisOrientation::CounterClockwise => "counter_clockwise",
    }
}
