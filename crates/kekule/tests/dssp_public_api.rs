use kekule::dssp::{
    self, DsspChainBreak, DsspError, DsspLimits, DsspOptions, DsspResource, DsspSecondaryStructure,
};
use kekule::mmcif::{self, MmcifInterpretOptions, MmcifModelSelection};
use kekule::topology::ResidueId;

const CRAMBIN_MMCIF: &str = include_str!("../../../benchmarks/corpora/smoke/data/rcsb/1CRN.cif");

fn crambin_model() -> kekule::structure::Model {
    let document = mmcif::parse_str(CRAMBIN_MMCIF).expect("checked-in RCSB 1CRN fixture parses");
    mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            model_selection: MmcifModelSelection::First,
            ..MmcifInterpretOptions::default()
        },
    )
    .expect("checked-in RCSB 1CRN fixture interprets")
    .into_model()
}

#[test]
fn dssp_matches_biopython_mkdssp_4_6_1_for_crambin() {
    let model = crambin_model();
    let result = dssp::assign(model.view(), DsspOptions::default())
        .expect("crambin has an analyzable protein backbone");
    let residues = result.residues().collect::<Vec<_>>();
    let codes = residues
        .iter()
        .map(|residue| residue.secondary_structure().code())
        .collect::<String>();

    assert_eq!(residues.len(), 46);
    let first_key: ResidueId = residues[0].key();
    assert!(result.get(first_key).is_some());
    assert_eq!(codes, " EE SSHHHHHHHHHHHTTT  HHHHHHHHS EE SSS   GGG  ");
    assert_eq!(residues[0].chain_break(), DsspChainBreak::NewChain);
    assert_eq!(
        residues[1].secondary_structure(),
        DsspSecondaryStructure::ExtendedStrand
    );
    assert!((residues[1].phi_degrees().expect("defined phi") - -107.8).abs() < 0.15);
    assert!((residues[1].psi_degrees().expect("defined psi") - 144.3).abs() < 0.15);

    let residue_three = residues[2];
    let strongest_acceptor = residue_three.acceptors()[0].expect("reference acceptor");
    let strongest_donor = residue_three.donors()[0].expect("reference donor");
    assert_eq!(
        result
            .get(strongest_acceptor.partner)
            .expect("acceptor partner is in the result")
            .source()
            .label_sequence_id,
        Some(33)
    );
    assert!((strongest_acceptor.energy_kcal_per_mol - -2.4).abs() < 0.05);
    assert_eq!(
        result
            .get(strongest_donor.partner)
            .expect("donor partner is in the result")
            .source()
            .label_sequence_id,
        Some(33)
    );
    assert!((strongest_donor.energy_kcal_per_mol - -2.8).abs() < 0.05);
}

#[test]
fn dssp_is_a_coordinate_snapshot_and_does_not_mutate_the_model() {
    let mut model = crambin_model();
    let before = model.clone();
    let assigned = dssp::assign(model.view(), DsspOptions::default()).expect("initial assignment");
    assert_eq!(model, before);

    let first_atom = model.topology().atom_ids()[0];
    let mut position = model.position(first_atom).expect("first position");
    position.x += 0.25;
    model
        .set_position(first_atom, position)
        .expect("finite coordinate update");

    assert_eq!(assigned.residues().count(), 46);
    assert_ne!(model, before);
}

#[test]
fn dssp_rejects_invalid_options_and_residue_limits() {
    let model = crambin_model();
    let options = DsspOptions {
        min_polyproline_stretch: 4,
        ..DsspOptions::default()
    };
    assert_eq!(
        dssp::assign(model.view(), options),
        Err(DsspError::InvalidPolyprolineStretch { value: 4 })
    );

    let options = DsspOptions {
        limits: DsspLimits {
            max_residues: 2,
            ..DsspLimits::default()
        },
        ..DsspOptions::default()
    };
    assert_eq!(
        dssp::assign(model.view(), options),
        Err(DsspError::ResourceLimitExceeded {
            resource: DsspResource::Residues,
            limit: 2,
        })
    );

    let options = DsspOptions {
        limits: DsspLimits {
            max_candidate_pairs: 0,
            ..DsspLimits::default()
        },
        ..DsspOptions::default()
    };
    assert_eq!(
        dssp::assign(model.view(), options),
        Err(DsspError::ResourceLimitExceeded {
            resource: DsspResource::CandidatePairs,
            limit: 0,
        })
    );

    let options = DsspOptions {
        limits: DsspLimits {
            max_ladders: 0,
            ..DsspLimits::default()
        },
        ..DsspOptions::default()
    };
    assert_eq!(
        dssp::assign(model.view(), options),
        Err(DsspError::ResourceLimitExceeded {
            resource: DsspResource::Ladders,
            limit: 0,
        })
    );
}

#[test]
fn dssp_rejects_coordinates_outside_its_spatial_index_range_without_panicking() {
    let mut model = crambin_model();
    let first_atom = model.topology().atom_ids()[0];
    let mut position = model.position(first_atom).expect("first position");
    position.x = f32::MAX as f64;
    model
        .set_position(first_atom, position)
        .expect("coordinate remains finite in the model");

    assert!(matches!(
        dssp::assign(model.view(), DsspOptions::default()),
        Err(DsspError::CoordinateOutOfRange {
            quantity: "backbone coordinate",
            ..
        })
    ));
}

#[test]
fn dssp_codes_cover_the_complete_dssp4_alphabet() {
    let codes = DsspSecondaryStructure::ALL.map(DsspSecondaryStructure::code);
    assert_eq!(codes, [' ', 'H', 'B', 'E', 'G', 'I', 'P', 'T', 'S']);
    for code in codes {
        assert_eq!(DsspSecondaryStructure::try_from(code).unwrap().code(), code);
    }
    assert_eq!(
        DsspSecondaryStructure::try_from('-').unwrap(),
        DsspSecondaryStructure::Loop
    );
    assert_eq!(
        DsspSecondaryStructure::try_from('X').unwrap_err().code(),
        'X'
    );
}

#[test]
fn dssp_analyzes_canonical_amino_acids_without_label_sequence_numbers() {
    let reference = dssp::assign(crambin_model().view(), DsspOptions::default()).unwrap();
    let source = CRAMBIN_MMCIF.replace("_atom_site.label_seq_id", "_atom_site.audit_sequence");
    let document = mmcif::parse_str(&source).unwrap();
    let model = mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            model_selection: MmcifModelSelection::First,
            ..MmcifInterpretOptions::default()
        },
    )
    .unwrap()
    .into_model();
    assert!(model
        .topology()
        .hierarchy()
        .residues()
        .all(|(_, residue)| residue.label_seq_id().is_none()));
    let assigned = dssp::assign(model.view(), DsspOptions::default()).unwrap();
    assert_eq!(assigned.statistics(), reference.statistics());
    for (actual, expected) in assigned.residues().zip(reference.residues()) {
        assert_eq!(
            actual.source().author_sequence_id,
            expected.source().author_sequence_id
        );
        assert_eq!(actual.secondary_structure(), expected.secondary_structure());
        assert_eq!(actual.phi_degrees(), expected.phi_degrees());
        assert_eq!(actual.psi_degrees(), expected.psi_degrees());
        assert_eq!(actual.kappa_degrees(), expected.kappa_degrees());
    }
}

#[test]
fn dssp_instance_membership_handles_many_reused_solvent_instances() {
    use kekule::core::AtomId;
    use kekule::geometry::Point3;
    use kekule::structure::Positions;
    use kekule::topology::{AtomSiteMetadata, InstanceAtomId};
    use kekule::units::{Quantity, ANGSTROM};

    let mut model = crambin_model().into_builder();
    let water = kekule::smiles::to_molecules("O").unwrap().pop().unwrap();
    let definition = model.add_molecule_definition(&water).unwrap();
    let water_position =
        Positions::new(Quantity::new(vec![Point3::new(100., 100., 100.)], ANGSTROM)).unwrap();
    let chain = model.hierarchy_mut().add_chain("water", None).unwrap();
    let mut ignored = Vec::new();
    let count = 10_000;
    for index in 0..count {
        let instance = model.add_instance(definition, &water_position).unwrap();
        if index % 11 == 0 {
            ignored.push(instance);
        } else {
            let residue = model
                .hierarchy_mut()
                .add_residue(chain, "HOH", None, None, None)
                .unwrap();
            model
                .hierarchy_mut()
                .add_atom_site(
                    residue,
                    InstanceAtomId::new(instance, AtomId::new(0)),
                    AtomSiteMetadata::default(),
                )
                .unwrap();
        }
    }
    let model = model.build().unwrap();
    let result = dssp::assign(model.view(), DsspOptions::default()).unwrap();
    assert_eq!(result.report().ignored_instances(), ignored);
    assert_eq!(result.statistics().analyzed_residues(), 46);
    assert_eq!(
        result.report().non_peptide_residues(),
        count - ignored.len()
    );
}

#[test]
fn dssp_numbering_fallback_is_independent_of_unlabeled_solvent_in_the_chain() {
    use kekule::core::AtomId;
    use kekule::geometry::Point3;
    use kekule::structure::{Model, Positions};
    use kekule::topology::{AtomSiteMetadata, InstanceAtomId, TopologyBuilder};
    use kekule::units::{Quantity, ANGSTROM};

    for (labels, authors) in [(true, false), (false, true), (false, false)] {
        let mut builder = TopologyBuilder::new();
        let molecule = kekule::smiles::to_molecules("NCC=O")
            .unwrap()
            .pop()
            .unwrap();
        let definition = builder.add_molecule_definition(&molecule).unwrap();
        let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
        let mut points = Vec::new();
        let mut keys = Vec::new();
        for sequence in [2, 1] {
            let instance = builder.add_instance(definition).unwrap();
            let residue = builder
                .hierarchy_mut()
                .add_residue(
                    chain,
                    "GLY",
                    labels.then_some(sequence),
                    authors.then(|| sequence.to_string()),
                    None,
                )
                .unwrap();
            keys.push(residue);
            for (atom, name) in molecule.atom_ids().zip(["N", "CA", "C", "O"]) {
                builder
                    .hierarchy_mut()
                    .add_atom_site(
                        residue,
                        InstanceAtomId::new(instance, atom),
                        AtomSiteMetadata {
                            label_atom_id: Some(name.to_owned()),
                            ..AtomSiteMetadata::default()
                        },
                    )
                    .unwrap();
            }
            let shift = f64::from(sequence - 1) * 3.;
            points.extend(
                [(0., 0.), (1.4, 0.), (2., 1.2), (3.1, 1.3)]
                    .map(|(x, y)| Point3::new(x + shift, y, 0.)),
            );
        }
        let water = kekule::smiles::to_molecules("O").unwrap().pop().unwrap();
        let instance = builder.add_molecule(&water).unwrap();
        let residue = builder
            .hierarchy_mut()
            .add_residue(chain, "HOH", None, None, None)
            .unwrap();
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, AtomId::new(0)),
                AtomSiteMetadata::default(),
            )
            .unwrap();
        points.push(Point3::new(100., 0., 0.));
        let model = Model::new(
            builder.build().unwrap(),
            Positions::new(Quantity::new(points, ANGSTROM)).unwrap(),
        )
        .unwrap();
        let result = dssp::assign(model.view(), DsspOptions::default()).unwrap();
        if labels || authors {
            keys.reverse();
        }
        assert_eq!(
            result
                .residues()
                .map(|residue| residue.key())
                .collect::<Vec<_>>(),
            keys
        );
        assert_eq!(result.report().non_peptide_residues(), 1);
    }
}
