use crate::core::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StereoPerceptionOptions {
    pub validate_existing: bool,
    pub detect_candidates: bool,
    pub assemble_source_marks: bool,
}

impl Default for StereoPerceptionOptions {
    fn default() -> Self {
        Self {
            validate_existing: true,
            detect_candidates: true,
            assemble_source_marks: true,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct StereoPerceptionReport {
    pub candidates: Vec<StereoCandidate>,
    pub issues: Vec<StereoPerceptionIssue>,
    pub assembled_elements: Vec<StereoElement>,
    pub created_elements: Vec<StereoElementId>,
}

impl StereoPerceptionReport {
    pub fn is_ok(&self) -> bool {
        self.issues.is_empty()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StereoCandidate {
    Tetrahedral {
        center: AtomId,
        carriers: Vec<StereoCarrier>,
    },
    DoubleBond {
        bond: BondId,
        left: AtomId,
        right: AtomId,
        left_carriers: Vec<StereoCarrier>,
        right_carriers: Vec<StereoCarrier>,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StereoPerceptionIssue {
    MissingStereoAtom {
        element: StereoElementId,
        atom: AtomId,
    },
    MissingStereoBond {
        element: StereoElementId,
        bond: BondId,
    },
    InvalidTetrahedralCarrierCount {
        element: StereoElementId,
        center: AtomId,
        carrier_count: usize,
    },
    DuplicateTetrahedralCarrier {
        element: StereoElementId,
        center: AtomId,
        carrier: StereoCarrier,
    },
    TetrahedralCarrierNotAdjacent {
        element: StereoElementId,
        center: AtomId,
        carrier: StereoCarrier,
    },
    TetrahedralHydrogenCarrierUnavailable {
        element: StereoElementId,
        center: AtomId,
    },
    InvalidDoubleBondOrder {
        element: StereoElementId,
        bond: BondId,
        order: BondOrder,
    },
    DoubleBondFocusMismatch {
        element: StereoElementId,
        bond: BondId,
        left: AtomId,
        right: AtomId,
    },
    DoubleBondCarrierIsFocusAtom {
        element: StereoElementId,
        endpoint: AtomId,
        carrier: AtomId,
    },
    DoubleBondCarrierNotAdjacent {
        element: StereoElementId,
        endpoint: AtomId,
        carrier: StereoCarrier,
    },
    DoubleBondHydrogenCarrierUnavailable {
        element: StereoElementId,
        endpoint: AtomId,
    },
    UnsupportedAxisElement {
        element: StereoElementId,
    },
    AmbiguousDirectionalBondMarks {
        double_bond: BondId,
        endpoint: AtomId,
        mark_count: usize,
    },
    UnpairedDirectionalBondMark {
        bond: BondId,
    },
    UnsupportedSourceBondMark {
        bond: BondId,
        kind: StereoBondMarkKind,
    },
    CouldNotCreateElement {
        message: String,
    },
}

pub fn validate_stereo(mol: &Molecule) -> StereoPerceptionReport {
    validate_stereo_with_options(mol, StereoPerceptionOptions::default())
}

pub fn validate_stereo_with_options(
    mol: &Molecule,
    options: StereoPerceptionOptions,
) -> StereoPerceptionReport {
    stereo_report(mol, options)
}

pub fn perceive_stereo(mol: &mut Molecule) -> StereoPerceptionReport {
    perceive_stereo_with_options(mol, StereoPerceptionOptions::default())
}

pub fn perceive_stereo_with_options(
    mol: &mut Molecule,
    options: StereoPerceptionOptions,
) -> StereoPerceptionReport {
    let mut report = stereo_report(mol, options);
    for element in report.assembled_elements.clone() {
        match mol.add_stereo_element(element) {
            Ok(id) => report.created_elements.push(id),
            Err(error) => report
                .issues
                .push(StereoPerceptionIssue::CouldNotCreateElement {
                    message: error.to_string(),
                }),
        }
    }
    mol.perception.stereo = ComputedState::Fresh;
    report
}

fn stereo_report(mol: &Molecule, options: StereoPerceptionOptions) -> StereoPerceptionReport {
    let mut report = StereoPerceptionReport::default();
    if options.validate_existing {
        validate_existing_elements(mol, &mut report.issues);
    }
    if options.detect_candidates {
        report.candidates.extend(tetrahedral_candidates(mol));
        report.candidates.extend(double_bond_candidates(mol));
    }
    if options.assemble_source_marks {
        report
            .assembled_elements
            .extend(assemble_directional_double_bonds(mol, &mut report.issues));
    }
    report
}

fn validate_existing_elements(mol: &Molecule, issues: &mut Vec<StereoPerceptionIssue>) {
    for (id, element) in mol.stereo_elements() {
        match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => validate_tetrahedral(mol, id, stereo, issues),
            StereoElementKind::DoubleBond(stereo) => validate_double_bond(mol, id, stereo, issues),
            StereoElementKind::Axis(_) => {
                issues.push(StereoPerceptionIssue::UnsupportedAxisElement { element: id })
            }
        }
    }
}

fn validate_tetrahedral(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &TetrahedralStereo,
    issues: &mut Vec<StereoPerceptionIssue>,
) {
    if mol.atom(stereo.center).is_err() {
        issues.push(StereoPerceptionIssue::MissingStereoAtom {
            element,
            atom: stereo.center,
        });
        return;
    }
    if stereo.carriers.len() != 4 {
        issues.push(StereoPerceptionIssue::InvalidTetrahedralCarrierCount {
            element,
            center: stereo.center,
            carrier_count: stereo.carriers.len(),
        });
    }
    let mut seen = Vec::<StereoCarrier>::new();
    for carrier in &stereo.carriers {
        if seen.contains(carrier) {
            issues.push(StereoPerceptionIssue::DuplicateTetrahedralCarrier {
                element,
                center: stereo.center,
                carrier: *carrier,
            });
        } else {
            seen.push(*carrier);
        }
        match carrier {
            StereoCarrier::Atom(atom) => {
                if mol.atom(*atom).is_err() {
                    issues.push(StereoPerceptionIssue::MissingStereoAtom {
                        element,
                        atom: *atom,
                    });
                } else if mol
                    .bond_between(stereo.center, *atom)
                    .ok()
                    .flatten()
                    .is_none()
                {
                    issues.push(StereoPerceptionIssue::TetrahedralCarrierNotAdjacent {
                        element,
                        center: stereo.center,
                        carrier: *carrier,
                    });
                }
            }
            StereoCarrier::ImplicitHydrogen => {
                if hydrogen_count(mol, stereo.center) == 0 {
                    issues.push(
                        StereoPerceptionIssue::TetrahedralHydrogenCarrierUnavailable {
                            element,
                            center: stereo.center,
                        },
                    );
                }
            }
        }
    }
}

fn validate_double_bond(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &DoubleBondStereo,
    issues: &mut Vec<StereoPerceptionIssue>,
) {
    let Ok(bond) = mol.bond(stereo.bond) else {
        issues.push(StereoPerceptionIssue::MissingStereoBond {
            element,
            bond: stereo.bond,
        });
        return;
    };
    if bond.order != BondOrder::Double {
        issues.push(StereoPerceptionIssue::InvalidDoubleBondOrder {
            element,
            bond: stereo.bond,
            order: bond.order,
        });
    }
    if !bond_connects(bond, stereo.left, stereo.right) {
        issues.push(StereoPerceptionIssue::DoubleBondFocusMismatch {
            element,
            bond: stereo.bond,
            left: stereo.left,
            right: stereo.right,
        });
    }
    validate_double_bond_carrier(
        mol,
        element,
        stereo.left,
        stereo.right,
        stereo.left_carrier,
        issues,
    );
    validate_double_bond_carrier(
        mol,
        element,
        stereo.right,
        stereo.left,
        stereo.right_carrier,
        issues,
    );
}

fn validate_double_bond_carrier(
    mol: &Molecule,
    element: StereoElementId,
    endpoint: AtomId,
    other_endpoint: AtomId,
    carrier: StereoCarrier,
    issues: &mut Vec<StereoPerceptionIssue>,
) {
    match carrier {
        StereoCarrier::Atom(atom) => {
            if atom == endpoint || atom == other_endpoint {
                issues.push(StereoPerceptionIssue::DoubleBondCarrierIsFocusAtom {
                    element,
                    endpoint,
                    carrier: atom,
                });
            } else if mol.atom(atom).is_err() {
                issues.push(StereoPerceptionIssue::MissingStereoAtom { element, atom });
            } else if mol.bond_between(endpoint, atom).ok().flatten().is_none() {
                issues.push(StereoPerceptionIssue::DoubleBondCarrierNotAdjacent {
                    element,
                    endpoint,
                    carrier,
                });
            }
        }
        StereoCarrier::ImplicitHydrogen => {
            if hydrogen_count(mol, endpoint) == 0 {
                issues.push(
                    StereoPerceptionIssue::DoubleBondHydrogenCarrierUnavailable {
                        element,
                        endpoint,
                    },
                );
            }
        }
    }
}

fn tetrahedral_candidates(mol: &Molecule) -> Vec<StereoCandidate> {
    let mut candidates = Vec::new();
    for (center, atom) in mol.atoms() {
        if atom.element.symbol() == "H" {
            continue;
        }
        let Ok(incident) = mol.incident_bonds(center) else {
            continue;
        };
        let mut atom_carriers = Vec::new();
        let mut single_bonded = true;
        for (_, bond) in incident {
            single_bonded &= bond.order == BondOrder::Single;
            atom_carriers.push(StereoCarrier::Atom(bond.other_atom(center)));
        }
        atom_carriers.sort_by_key(carrier_key);
        let hydrogens = hydrogen_count(mol, center);
        if single_bonded && hydrogens <= 1 && atom_carriers.len() + usize::from(hydrogens) == 4 {
            if hydrogens == 1 {
                atom_carriers.push(StereoCarrier::ImplicitHydrogen);
            }
            candidates.push(StereoCandidate::Tetrahedral {
                center,
                carriers: atom_carriers,
            });
        }
    }
    candidates
}

fn double_bond_candidates(mol: &Molecule) -> Vec<StereoCandidate> {
    let mut candidates = Vec::new();
    for (bond_id, bond) in mol.bonds() {
        if bond.order != BondOrder::Double {
            continue;
        }
        let left = bond.a();
        let right = bond.b();
        let left_carriers = double_bond_endpoint_carriers(mol, left, right, bond_id);
        let right_carriers = double_bond_endpoint_carriers(mol, right, left, bond_id);
        if !left_carriers.is_empty() && !right_carriers.is_empty() {
            candidates.push(StereoCandidate::DoubleBond {
                bond: bond_id,
                left,
                right,
                left_carriers,
                right_carriers,
            });
        }
    }
    candidates
}

fn double_bond_endpoint_carriers(
    mol: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    focus_bond: BondId,
) -> Vec<StereoCarrier> {
    let mut carriers = Vec::new();
    if let Ok(incident) = mol.incident_bonds(endpoint) {
        for (bond_id, bond) in incident {
            if bond_id == focus_bond || bond.order != BondOrder::Single {
                continue;
            }
            let other = bond.other_atom(endpoint);
            if other != other_endpoint {
                carriers.push(StereoCarrier::Atom(other));
            }
        }
    }
    carriers.sort_by_key(carrier_key);
    if hydrogen_count(mol, endpoint) == 1 {
        carriers.push(StereoCarrier::ImplicitHydrogen);
    }
    carriers
}

fn assemble_directional_double_bonds(
    mol: &Molecule,
    issues: &mut Vec<StereoPerceptionIssue>,
) -> Vec<StereoElement> {
    let mut assembled = Vec::new();
    let mut used_marks = Vec::<BondId>::new();
    for (bond_id, bond) in mol.bonds() {
        if bond.order != BondOrder::Double || has_double_bond_element(mol, bond_id) {
            continue;
        }
        let left = bond.a();
        let right = bond.b();
        let left_marks = directional_marks_for_endpoint(mol, left, bond_id);
        let right_marks = directional_marks_for_endpoint(mol, right, bond_id);
        if left_marks.len() > 1 {
            issues.push(StereoPerceptionIssue::AmbiguousDirectionalBondMarks {
                double_bond: bond_id,
                endpoint: left,
                mark_count: left_marks.len(),
            });
        }
        if right_marks.len() > 1 {
            issues.push(StereoPerceptionIssue::AmbiguousDirectionalBondMarks {
                double_bond: bond_id,
                endpoint: right,
                mark_count: right_marks.len(),
            });
        }
        let ([left_mark], [right_mark]) = (left_marks.as_slice(), right_marks.as_slice()) else {
            continue;
        };
        let orientation = if left_mark.mark.kind == right_mark.mark.kind {
            DoubleBondOrientation::Together
        } else {
            DoubleBondOrientation::Opposite
        };
        used_marks.push(left_mark.bond);
        used_marks.push(right_mark.bond);
        assembled.push(StereoElement::specified(
            StereoElementKind::DoubleBond(DoubleBondStereo {
                bond: bond_id,
                left,
                right,
                left_carrier: StereoCarrier::Atom(left_mark.carrier),
                right_carrier: StereoCarrier::Atom(right_mark.carrier),
                orientation,
            }),
            common_source(left_mark.mark.source, right_mark.mark.source),
        ));
    }

    for mark in mol.stereo_bond_marks() {
        match mark.kind {
            StereoBondMarkKind::DirectionalUp | StereoBondMarkKind::DirectionalDown => {
                if !used_marks.contains(&mark.bond) {
                    issues.push(StereoPerceptionIssue::UnpairedDirectionalBondMark {
                        bond: mark.bond,
                    });
                }
            }
            StereoBondMarkKind::WedgeUp
            | StereoBondMarkKind::WedgeDown
            | StereoBondMarkKind::WedgeEither
            | StereoBondMarkKind::DoubleBondEither => {
                issues.push(StereoPerceptionIssue::UnsupportedSourceBondMark {
                    bond: mark.bond,
                    kind: mark.kind,
                });
            }
        }
    }
    assembled
}

#[derive(Clone, Copy)]
struct EndpointMark<'a> {
    bond: BondId,
    carrier: AtomId,
    mark: &'a StereoBondMark,
}

fn directional_marks_for_endpoint(
    mol: &Molecule,
    endpoint: AtomId,
    focus_bond: BondId,
) -> Vec<EndpointMark<'_>> {
    let mut marks = Vec::new();
    let Ok(incident) = mol.incident_bonds(endpoint) else {
        return marks;
    };
    for (bond_id, bond) in incident {
        if bond_id == focus_bond || bond.order != BondOrder::Single {
            continue;
        }
        let Some(mark) = mol.stereo_bond_mark(bond_id) else {
            continue;
        };
        if matches!(
            mark.kind,
            StereoBondMarkKind::DirectionalUp | StereoBondMarkKind::DirectionalDown
        ) {
            marks.push(EndpointMark {
                bond: bond_id,
                carrier: bond.other_atom(endpoint),
                mark,
            });
        }
    }
    marks.sort_by_key(|mark| (mark.bond, mark.carrier));
    marks
}

fn has_double_bond_element(mol: &Molecule, bond: BondId) -> bool {
    mol.stereo_elements().any(|(_, element)| {
        matches!(
            &element.kind,
            StereoElementKind::DoubleBond(stereo) if stereo.bond == bond
        )
    })
}

fn hydrogen_count(mol: &Molecule, atom: AtomId) -> u8 {
    let Ok(atom) = mol.atom(atom) else {
        return 0;
    };
    atom.explicit_hydrogens
        .saturating_add(atom.implicit_hydrogens.unwrap_or(0))
}

fn carrier_key(carrier: &StereoCarrier) -> (u8, u32) {
    match carrier {
        StereoCarrier::Atom(atom) => (0, atom.raw()),
        StereoCarrier::ImplicitHydrogen => (1, u32::MAX),
    }
}

fn bond_connects(bond: &Bond, a: AtomId, b: AtomId) -> bool {
    (bond.a() == a && bond.b() == b) || (bond.a() == b && bond.b() == a)
}

fn common_source(left: StereoSource, right: StereoSource) -> StereoSource {
    if left == right {
        left
    } else {
        StereoSource::User
    }
}
