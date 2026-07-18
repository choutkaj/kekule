use molecules::dssp::{
    self, DsspChainBreak, DsspError, DsspLimits, DsspOptions, DsspResource, DsspSecondaryStructure,
};
use molecules::mmcif::{self, MmcifInterpretOptions, MmcifModelSelection, MmcifParseOptions};

const CRAMBIN_MMCIF: &str = include_str!("../../../validation/corpora/smoke/data/rcsb/1CRN.cif");

fn crambin_model() -> molecules::modeling::Model {
    let document = mmcif::parse_str(CRAMBIN_MMCIF, MmcifParseOptions::default())
        .expect("checked-in RCSB 1CRN fixture parses");
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
    let result = dssp::assign(&crambin_model(), DsspOptions::default())
        .expect("crambin has an analyzable protein backbone");
    let residues = result.residues().collect::<Vec<_>>();
    let codes = residues
        .iter()
        .map(|residue| residue.secondary_structure().code())
        .collect::<String>();

    assert_eq!(residues.len(), 46);
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
    assert_eq!(strongest_acceptor.partner.residue().raw(), 32);
    assert!((strongest_acceptor.energy_kcal_per_mol - -2.4).abs() < 0.05);
    assert_eq!(strongest_donor.partner.residue().raw(), 32);
    assert!((strongest_donor.energy_kcal_per_mol - -2.8).abs() < 0.05);
}

#[test]
fn dssp_is_a_coordinate_snapshot_and_does_not_mutate_the_model() {
    let mut model = crambin_model();
    let before = model.clone();
    let assigned = dssp::assign(&model, DsspOptions::default()).expect("initial assignment");
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
        dssp::assign(&model, options),
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
        dssp::assign(&model, options),
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
        dssp::assign(&model, options),
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
        dssp::assign(&model, options),
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
        dssp::assign(&model, DsspOptions::default()),
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
