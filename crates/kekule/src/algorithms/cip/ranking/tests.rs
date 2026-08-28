use super::*;

fn node_priority(atomic_number: u8, rule1b: u32, isotope: u16) -> NodePriority {
    NodePriority {
        atomic_number: AtomicNumberFraction::element(atomic_number),
        rule1b,
        rule2_mass: if isotope == 0 {
            Rule2Mass::natural(atomic_number)
        } else {
            Rule2Mass::isotope(isotope)
        },
        descriptor: None,
        rule6_atom: None,
    }
}

fn node_priority_with_descriptor(descriptor: StereoDescriptor) -> NodePriority {
    NodePriority {
        descriptor: Some(descriptor),
        ..node_priority(6, 0, 0)
    }
}

fn node_priority_with_rule6_atom(atomic_number: u8, atom: AtomId) -> NodePriority {
    NodePriority {
        rule6_atom: Some(atom),
        ..node_priority(atomic_number, 0, 0)
    }
}

fn one_node_signature(rule1b: u32, isotope: u16) -> LigandTree {
    LigandTree {
        priority: node_priority(6, rule1b, isotope),
        children: Vec::new(),
    }
}

fn signature(root: LigandTree) -> LigandSignature {
    LigandSignature { root }
}

#[test]
fn rule1b_ring_duplicate_priority_is_applied_before_isotope_priority() {
    let ring_duplicate = signature(one_node_signature(u32::MAX, 0));
    let isotope = signature(one_node_signature(0, 13));

    assert_eq!(ring_duplicate.compare(&isotope), Ordering::Greater);
    assert_eq!(isotope.compare(&ring_duplicate), Ordering::Less);
}

#[test]
fn rule2_compares_indicated_isotopes_against_natural_atomic_weight() {
    let natural_hydrogen = node_priority(1, 0, 0);
    let protium = node_priority(1, 0, 1);
    let deuterium = node_priority(1, 0, 2);

    assert_eq!(
        natural_hydrogen.compare_by_rule(&protium, SequenceRule::Rule2, None),
        Ordering::Greater
    );
    assert_eq!(
        deuterium.compare_by_rule(&natural_hydrogen, SequenceRule::Rule2, None),
        Ordering::Greater
    );

    let natural_carbon = node_priority(6, 0, 0);
    let carbon_12 = node_priority(6, 0, 12);
    let carbon_13 = node_priority(6, 0, 13);

    assert_eq!(
        natural_carbon.compare_by_rule(&carbon_12, SequenceRule::Rule2, None),
        Ordering::Greater
    );
    assert_eq!(
        carbon_13.compare_by_rule(&natural_carbon, SequenceRule::Rule2, None),
        Ordering::Greater
    );

    let another_natural_hydrogen = node_priority(1, 0, 0);
    assert_eq!(
        natural_hydrogen.compare_by_rule(&another_natural_hydrogen, SequenceRule::Rule2, None),
        Ordering::Equal
    );
}

#[test]
fn rule1b_prefers_ring_duplicate_whose_reference_is_closer_to_root() {
    let root_reference = signature(one_node_signature(u32::MAX, 0));
    let deeper_reference = signature(one_node_signature(u32::MAX - 2, 0));

    assert_eq!(root_reference.compare(&deeper_reference), Ordering::Greater);
    assert_eq!(deeper_reference.compare(&root_reference), Ordering::Less);
}

#[test]
fn rule4a_prefers_uppercase_sequence_descriptors_over_pseudo_descriptors() {
    let uppercase = signature(LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::R),
        children: Vec::new(),
    });
    let uppercase_axis = signature(LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::M),
        children: Vec::new(),
    });
    let sequence = signature(LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::SeqCis),
        children: Vec::new(),
    });
    let pseudo = signature(LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::LowerR),
        children: Vec::new(),
    });
    let pseudo_axis = signature(LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::LowerM),
        children: Vec::new(),
    });
    let unlabeled = signature(one_node_signature(0, 0));

    assert_eq!(uppercase.compare(&pseudo), Ordering::Greater);
    assert_eq!(uppercase_axis.compare(&pseudo_axis), Ordering::Greater);
    assert_eq!(sequence.compare(&pseudo), Ordering::Greater);
    assert_eq!(sequence.compare(&uppercase), Ordering::Equal);
    assert_eq!(pseudo_axis.compare(&pseudo), Ordering::Equal);
    assert_eq!(pseudo.compare(&unlabeled), Ordering::Greater);
    assert_eq!(unlabeled.compare(&uppercase), Ordering::Less);
}

#[test]
fn rule4b_reference_descriptors_use_first_equivalent_descriptor_level() {
    let majority_r = LigandTree {
        priority: node_priority(6, 0, 0),
        children: vec![
            LigandTree {
                priority: node_priority_with_descriptor(StereoDescriptor::R),
                children: Vec::new(),
            },
            LigandTree {
                priority: node_priority_with_descriptor(StereoDescriptor::M),
                children: Vec::new(),
            },
            LigandTree {
                priority: node_priority_with_descriptor(StereoDescriptor::S),
                children: Vec::new(),
            },
        ],
    };
    let tied = LigandTree {
        priority: node_priority(6, 0, 0),
        children: vec![
            LigandTree {
                priority: node_priority_with_descriptor(StereoDescriptor::R),
                children: Vec::new(),
            },
            LigandTree {
                priority: node_priority_with_descriptor(StereoDescriptor::P),
                children: Vec::new(),
            },
        ],
    };

    assert_eq!(
        majority_r.rule4b_reference_descriptors(),
        vec![DescriptorRef::R]
    );
    assert_eq!(
        tied.rule4b_reference_descriptors(),
        vec![DescriptorRef::R, DescriptorRef::S]
    );
}

#[test]
fn rule4b_fixed_reference_prefers_like_descriptor_families() {
    let r_ligand = LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::R),
        children: Vec::new(),
    };
    let s_ligand = LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::S),
        children: Vec::new(),
    };

    assert_eq!(
        r_ligand.compare_with_reference(&s_ligand, DescriptorRef::R),
        Ordering::Greater
    );
    assert_eq!(
        r_ligand.compare_with_reference(&s_ligand, DescriptorRef::S),
        Ordering::Less
    );

    let seq_cis = LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::SeqCis),
        children: Vec::new(),
    };
    let seq_trans = LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::SeqTrans),
        children: Vec::new(),
    };

    assert_eq!(
        seq_cis.compare_with_reference(&seq_trans, DescriptorRef::R),
        Ordering::Greater
    );
    assert_eq!(
        seq_cis.compare_with_reference(&seq_trans, DescriptorRef::S),
        Ordering::Less
    );
}

#[test]
fn rule4c_prefers_lower_r_and_lower_m_over_lower_s_and_lower_p() {
    let lower_r = signature(LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::LowerR),
        children: Vec::new(),
    });
    let lower_s = signature(LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::LowerS),
        children: Vec::new(),
    });
    let lower_m = signature(LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::LowerM),
        children: Vec::new(),
    });
    let lower_p = signature(LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::LowerP),
        children: Vec::new(),
    });

    assert_eq!(lower_r.compare(&lower_s), Ordering::Greater);
    assert_eq!(lower_s.compare(&lower_r), Ordering::Less);
    assert_eq!(lower_m.compare(&lower_p), Ordering::Greater);
    assert_eq!(lower_p.compare(&lower_m), Ordering::Less);
    assert_eq!(lower_m.compare(&lower_r), Ordering::Equal);
    assert_eq!(lower_p.compare(&lower_s), Ordering::Equal);
}

#[test]
fn rule5_descriptor_pairing_prefers_like_pairs_and_marks_pseudo_asymmetric_ordering() {
    let r_ligand = signature(LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::R),
        children: Vec::new(),
    });
    let s_ligand = signature(LigandTree {
        priority: node_priority_with_descriptor(StereoDescriptor::S),
        children: Vec::new(),
    });

    let comparison = r_ligand.compare_with_flags(&s_ligand);

    assert_eq!(comparison.ordering, Ordering::Greater);
    assert!(comparison.pseudo_asymmetric);

    let comparison = s_ligand.compare_with_flags(&r_ligand);

    assert_eq!(comparison.ordering, Ordering::Less);
    assert!(comparison.pseudo_asymmetric);
}

#[test]
fn rule6_prefers_nodes_matching_the_selected_reference_atom() {
    let reference = AtomId::new(7);
    let reference_ligand = signature(LigandTree {
        priority: node_priority_with_rule6_atom(6, reference),
        children: Vec::new(),
    });
    let other_ligand = signature(LigandTree {
        priority: node_priority_with_rule6_atom(6, AtomId::new(8)),
        children: Vec::new(),
    });

    let comparison = reference_ligand.compare_with_rule6_reference(&other_ligand, Some(reference));

    assert_eq!(comparison.ordering, Ordering::Greater);
    assert!(!comparison.pseudo_asymmetric);
    assert_eq!(reference_ligand.compare(&other_ligand), Ordering::Equal);
}

#[test]
fn rule6_tetrahedral_retry_resolves_two_equivalent_partitions() {
    let carrier_a = AtomId::new(0);
    let carrier_b = AtomId::new(1);
    let carrier_c = AtomId::new(2);
    let carrier_d = AtomId::new(3);
    let reference_child = LigandTree {
        priority: node_priority_with_rule6_atom(1, carrier_b),
        children: Vec::new(),
    };
    let other_child = LigandTree {
        priority: node_priority_with_rule6_atom(1, AtomId::new(99)),
        children: Vec::new(),
    };
    let signatures = vec![
        (
            StereoCarrier::Atom(carrier_a),
            signature(LigandTree {
                priority: node_priority_with_rule6_atom(8, carrier_a),
                children: Vec::new(),
            }),
        ),
        (
            StereoCarrier::Atom(carrier_b),
            signature(LigandTree {
                priority: node_priority_with_rule6_atom(8, carrier_b),
                children: Vec::new(),
            }),
        ),
        (
            StereoCarrier::Atom(carrier_c),
            signature(LigandTree {
                priority: node_priority_with_rule6_atom(6, carrier_c),
                children: vec![reference_child],
            }),
        ),
        (
            StereoCarrier::Atom(carrier_d),
            signature(LigandTree {
                priority: node_priority_with_rule6_atom(6, carrier_d),
                children: vec![other_child],
            }),
        ),
    ];

    let ranked = rank_tetrahedral_signatures_with_rule6(
        &crate::core::MoleculeEditor::new(),
        StereoElementId::new(0),
        AtomId::new(0),
        &signatures,
        TetrahedralOrientation::Clockwise,
        false,
    )
    .expect("Rule 6 should resolve paired partitions");

    assert_eq!(
        ranked.carriers,
        vec![
            StereoCarrier::Atom(carrier_b),
            StereoCarrier::Atom(carrier_a),
            StereoCarrier::Atom(carrier_c),
            StereoCarrier::Atom(carrier_d),
        ]
    );
    assert!(!ranked.pseudo_asymmetric_ordering);
}

#[test]
fn rule6_tetrahedral_retry_rejects_parity_unstable_two_partition_rankings() {
    let carrier_a = AtomId::new(0);
    let carrier_b = AtomId::new(1);
    let carrier_c = AtomId::new(2);
    let carrier_d = AtomId::new(3);
    let signatures = vec![
        (
            StereoCarrier::Atom(carrier_a),
            signature(LigandTree {
                priority: node_priority_with_rule6_atom(8, carrier_a),
                children: Vec::new(),
            }),
        ),
        (
            StereoCarrier::Atom(carrier_b),
            signature(LigandTree {
                priority: node_priority_with_rule6_atom(8, carrier_b),
                children: Vec::new(),
            }),
        ),
        (
            StereoCarrier::Atom(carrier_c),
            signature(LigandTree {
                priority: node_priority_with_rule6_atom(7, carrier_c),
                children: Vec::new(),
            }),
        ),
        (
            StereoCarrier::Atom(carrier_d),
            signature(LigandTree {
                priority: node_priority_with_rule6_atom(6, carrier_d),
                children: Vec::new(),
            }),
        ),
    ];
    let element = StereoElementId::new(0);

    let issue = rank_tetrahedral_signatures_with_rule6(
        &crate::core::MoleculeEditor::new(),
        element,
        AtomId::new(0),
        &signatures,
        TetrahedralOrientation::Clockwise,
        false,
    )
    .expect_err("odd successful reference permutations must remain unresolved");

    assert_eq!(issue, CipAssignmentIssue::UnresolvedPriority { element });
}

fn s4_rule6_signatures(
    child_reference_counts: [[usize; 4]; 4],
) -> Vec<(StereoCarrier, LigandSignature)> {
    let carriers = [
        AtomId::new(0),
        AtomId::new(1),
        AtomId::new(2),
        AtomId::new(3),
    ];
    carriers
        .iter()
        .copied()
        .enumerate()
        .map(|(carrier_index, carrier)| {
            let mut children = Vec::new();
            for (reference_index, reference) in carriers.iter().copied().enumerate() {
                for _ in 0..child_reference_counts[carrier_index][reference_index] {
                    children.push(LigandTree {
                        priority: node_priority_with_rule6_atom(1, reference),
                        children: Vec::new(),
                    });
                }
            }
            (
                StereoCarrier::Atom(carrier),
                signature(LigandTree {
                    priority: node_priority_with_rule6_atom(6, carrier),
                    children,
                }),
            )
        })
        .collect()
}

#[test]
fn rule6_s4_retry_accepts_parity_stable_reference_rankings() {
    let signatures = s4_rule6_signatures([[0, 2, 0, 2], [2, 0, 2, 0], [1, 0, 2, 1], [0, 1, 1, 2]]);

    let ranked = rank_tetrahedral_signatures_with_rule6(
        &crate::core::MoleculeEditor::new(),
        StereoElementId::new(0),
        AtomId::new(0),
        &signatures,
        TetrahedralOrientation::Clockwise,
        false,
    )
    .expect("Rule 6 should accept parity-stable S4 rankings");

    assert_eq!(
        ranked.carriers,
        vec![
            StereoCarrier::Atom(AtomId::new(0)),
            StereoCarrier::Atom(AtomId::new(1)),
            StereoCarrier::Atom(AtomId::new(2)),
            StereoCarrier::Atom(AtomId::new(3)),
        ]
    );
    assert!(!ranked.pseudo_asymmetric_ordering);
}

#[test]
fn rule6_s4_retry_rejects_parity_unstable_reference_rankings() {
    let element = StereoElementId::new(0);
    let signatures = s4_rule6_signatures([[0, 2, 1, 1], [2, 0, 0, 2], [1, 1, 2, 0], [0, 0, 2, 2]]);

    let issue = rank_tetrahedral_signatures_with_rule6(
        &crate::core::MoleculeEditor::new(),
        element,
        AtomId::new(0),
        &signatures,
        TetrahedralOrientation::Clockwise,
        false,
    )
    .expect_err("odd reference permutations must remain unresolved");

    assert_eq!(issue, CipAssignmentIssue::UnresolvedPriority { element });
}

#[test]
fn ligand_tree_compares_highest_priority_branch_before_lower_siblings() {
    let oxygen_to_carbon = LigandTree {
        priority: node_priority(8, 0, 0),
        children: vec![one_node_signature(0, 0)],
    };
    let oxygen_to_hydrogen = LigandTree {
        priority: node_priority(8, 0, 0),
        children: vec![LigandTree {
            priority: node_priority(1, 0, 0),
            children: Vec::new(),
        }],
    };
    let carbon_to_nitrogen = LigandTree {
        priority: node_priority(6, 0, 0),
        children: vec![LigandTree {
            priority: node_priority(7, 0, 0),
            children: Vec::new(),
        }],
    };
    let carbon_to_oxygen = LigandTree {
        priority: node_priority(6, 0, 0),
        children: vec![LigandTree {
            priority: node_priority(8, 0, 0),
            children: Vec::new(),
        }],
    };

    let left = signature(LigandTree {
        priority: node_priority(6, 0, 0),
        children: vec![oxygen_to_carbon, carbon_to_oxygen],
    });
    let right = signature(LigandTree {
        priority: node_priority(6, 0, 0),
        children: vec![oxygen_to_hydrogen, carbon_to_nitrogen],
    });

    assert_eq!(left.compare(&right), Ordering::Greater);
}

#[test]
fn ligand_tree_compares_immediate_sibling_list_before_recursing() {
    let oxygen_to_hydrogen = LigandTree {
        priority: node_priority(8, 0, 0),
        children: vec![LigandTree {
            priority: node_priority(1, 0, 0),
            children: Vec::new(),
        }],
    };
    let oxygen_to_phosphorus = LigandTree {
        priority: node_priority(8, 0, 0),
        children: vec![LigandTree {
            priority: node_priority(15, 0, 0),
            children: Vec::new(),
        }],
    };
    let carbon = one_node_signature(0, 0);
    let hydrogen = LigandTree {
        priority: node_priority(1, 0, 0),
        children: Vec::new(),
    };

    let left = signature(LigandTree {
        priority: node_priority(6, 0, 0),
        children: vec![oxygen_to_hydrogen, carbon],
    });
    let right = signature(LigandTree {
        priority: node_priority(6, 0, 0),
        children: vec![oxygen_to_phosphorus, hydrogen],
    });

    assert_eq!(left.compare(&right), Ordering::Greater);
}

#[test]
fn duplicate_nodes_have_no_isotope_priority() {
    let mut mol = crate::core::MoleculeEditor::new();
    let mut isotope = Atom::new(Element::from_symbol("C").expect("carbon"));
    isotope.isotope = Some(13);
    let atom = mol.add_atom(isotope).expect("atom identifier capacity");

    let normal = LigandNode::Atom {
        atom,
        previous: None,
        path: vec![atom],
        duplicate: None,
        terminal: false,
    };
    let duplicate = LigandNode::Atom {
        atom,
        previous: None,
        path: Vec::new(),
        duplicate: Some(DuplicateNode::Bond {
            atomic_number: None,
        }),
        terminal: true,
    };

    assert_eq!(normal.rule2_mass(&mol), Rule2Mass::isotope(13));
    assert_eq!(duplicate.rule2_mass(&mol), Rule2Mass::ZERO);
}

#[test]
fn rule1a_uses_mancude_fractional_atomic_numbers_for_bond_duplicates() {
    let mut mol = crate::core::MoleculeEditor::new();
    let atoms = (0..6)
        .map(|index| {
            let symbol = if index == 3 { "N" } else { "C" };
            mol.add_atom(Atom::new(Element::from_symbol(symbol).expect("element")))
                .expect("atom identifier capacity")
        })
        .collect::<Vec<_>>();
    for (left, right, order) in [
        (0, 1, BondOrder::Double),
        (1, 2, BondOrder::Single),
        (2, 3, BondOrder::Double),
        (3, 4, BondOrder::Single),
        (4, 5, BondOrder::Double),
        (5, 0, BondOrder::Single),
    ] {
        mol.add_bond(atoms[left], atoms[right], order)
            .expect("ring bond");
    }
    for (index, atom) in atoms.iter().copied().enumerate() {
        mol.set_implicit_hydrogens(atom, if index == 3 { 0 } else { 1 });
    }

    let cip_bond_orders = CipBondOrders::new(&mol, false);
    let fractions = cip_atomic_number_fractions(&mol, &cip_bond_orders);

    assert_eq!(
        fractions[atoms[2].index()],
        AtomicNumberFraction::new(13, 2)
    );
    assert_eq!(fractions[atoms[3].index()], AtomicNumberFraction::new(6, 1));
    assert_eq!(
        fractions[atoms[4].index()],
        AtomicNumberFraction::new(13, 2)
    );
    assert_eq!(fractions[atoms[0].index()], AtomicNumberFraction::new(6, 1));

    let node = LigandNode::Atom {
        atom: atoms[2],
        previous: Some(atoms[1]),
        path: vec![atoms[1], atoms[2]],
        duplicate: None,
        terminal: false,
    };
    let mut next = Vec::new();
    node.extend(&mol, &fractions, &cip_bond_orders, &mut next);

    let normal_nitrogen = next
        .iter()
        .find(|child| {
            matches!(
                child,
                LigandNode::Atom {
                    atom,
                    duplicate: None,
                    ..
                } if *atom == atoms[3]
            )
        })
        .expect("normal nitrogen child");
    let duplicate_nitrogen = next
        .iter()
        .find(|child| {
            matches!(
                child,
                LigandNode::Atom {
                    atom,
                    duplicate: Some(DuplicateNode::Bond { .. }),
                    ..
                } if *atom == atoms[3]
            )
        })
        .expect("duplicate nitrogen child");
    let element = StereoElementId::new(0);
    let descriptor_context = DescriptorContext::new(element, AuxiliaryDescriptorMode::Disabled);
    let build_context = LigandBuildContext {
        mol: &mol,
        element,
        descriptor_context: &descriptor_context,
        options: CipAssignmentOptions::default(),
        atomic_number_fractions: &fractions,
        cip_bond_orders: &cip_bond_orders,
    };

    assert_eq!(
        normal_nitrogen.priority(&build_context).atomic_number,
        AtomicNumberFraction::element(7)
    );
    assert_eq!(
        duplicate_nitrogen.priority(&build_context).atomic_number,
        AtomicNumberFraction::new(13, 2)
    );
}

#[test]
fn higher_order_bond_expansion_creates_terminal_duplicate_nodes() {
    let mut mol = crate::core::MoleculeEditor::new();
    let root = mol
        .add_atom(Atom::new(Element::from_symbol("C").expect("carbon")))
        .expect("atom identifier capacity");
    let carbon = mol
        .add_atom(Atom::new(Element::from_symbol("C").expect("carbon")))
        .expect("atom identifier capacity");
    let oxygen = mol
        .add_atom(Atom::new(Element::from_symbol("O").expect("oxygen")))
        .expect("atom identifier capacity");
    mol.add_bond(root, carbon, BondOrder::Single)
        .expect("root bond");
    mol.add_bond(carbon, oxygen, BondOrder::Double)
        .expect("double bond");

    let node = LigandNode::Atom {
        atom: carbon,
        previous: Some(root),
        path: vec![root, carbon],
        duplicate: None,
        terminal: false,
    };
    let mut next = Vec::new();
    let cip_bond_orders = CipBondOrders::new(&mol, false);
    let fractions = cip_atomic_number_fractions(&mol, &cip_bond_orders);
    node.extend(&mol, &fractions, &cip_bond_orders, &mut next);

    assert_eq!(next.len(), 2);
    assert!(next.contains(&LigandNode::Atom {
        atom: oxygen,
        previous: Some(carbon),
        path: vec![root, carbon, oxygen],
        duplicate: None,
        terminal: false,
    }));
    assert!(next.contains(&LigandNode::Atom {
        atom: oxygen,
        previous: Some(carbon),
        path: Vec::new(),
        duplicate: Some(DuplicateNode::Bond {
            atomic_number: None,
        }),
        terminal: true,
    }));
}

#[test]
fn negative_fractional_atoms_create_duplicate_nodes() {
    let mut mol = crate::core::MoleculeEditor::new();
    let atoms = (0..5)
        .map(|index| {
            let mut atom = Atom::new(Element::from_symbol("C").expect("carbon"));
            if index == 2 {
                atom.formal_charge = -1;
            }
            mol.add_atom(atom).expect("atom identifier capacity")
        })
        .collect::<Vec<_>>();
    for (left, right, order) in [
        (0, 1, BondOrder::Double),
        (1, 2, BondOrder::Single),
        (2, 3, BondOrder::Single),
        (3, 4, BondOrder::Double),
        (4, 0, BondOrder::Single),
    ] {
        mol.add_bond(atoms[left], atoms[right], order)
            .expect("ring bond");
    }
    for atom in &atoms {
        mol.set_implicit_hydrogens(*atom, 1);
    }

    let mut fractions = vec![AtomicNumberFraction::element(6); mol.graph.atom_slot_count()];
    fractions[atoms[2].index()] = AtomicNumberFraction::new(13, 2);
    let cip_bond_orders = CipBondOrders::new(&mol, false);

    let node = LigandNode::Atom {
        atom: atoms[2],
        previous: Some(atoms[1]),
        path: vec![atoms[1], atoms[2]],
        duplicate: None,
        terminal: false,
    };
    let mut next = Vec::new();
    node.extend(&mol, &fractions, &cip_bond_orders, &mut next);

    assert!(next.iter().any(|child| matches!(
        child,
        LigandNode::Atom {
            atom,
            duplicate: Some(DuplicateNode::Bond { .. }),
            terminal: true,
            ..
        } if *atom == atoms[3]
    )));
}
