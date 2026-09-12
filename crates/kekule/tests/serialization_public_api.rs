use std::sync::Arc;
use std::{io, io::Write};

use kekule::geometry::Point3;
use kekule::mmcif::{
    self, MmcifEnsembleInterpretOptions, MmcifEntityClassifications, MmcifEntityKind,
    MmcifWriteOptions,
};
use kekule::molfile::{self, MolWriteErrorKind, MolfileWriteOptions, MolfileWriteVersion};
use kekule::sdf::{self, SdfDataField, SdfRecordInterpretation, SdfWriteOptions};
use kekule::smiles::{self, SmilesWriteMode, SmilesWriteOptions};
use kekule::structure::{Ensemble, Model, Positions};
use kekule::topology::{Topology, TopologyBuilder};
use kekule::units::{Quantity, ANGSTROM};

fn molecule(source: &str) -> kekule::core::Molecule {
    let mut molecules = smiles::to_molecules(source).expect("SMILES interprets");
    assert_eq!(molecules.len(), 1);
    molecules.pop().unwrap()
}

fn model(source: &str, points: &[[f64; 3]]) -> Model {
    let topology = smiles::to_topology(source).expect("topology interprets");
    let positions = Positions::new(Quantity::new(
        points
            .iter()
            .map(|point| Point3::new(point[0], point[1], point[2]))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))
    .unwrap();
    Model::new(topology, positions).unwrap()
}

fn classifications(model: &Model) -> MmcifEntityClassifications {
    let mut classifications = MmcifEntityClassifications::new();
    for (instance, _) in model.topology().instances() {
        classifications
            .insert(instance, MmcifEntityKind::NonPolymer)
            .unwrap();
    }
    classifications
}

struct BoundedChunkWriter {
    bytes: Vec<u8>,
    maximum_chunk: usize,
}

impl Write for BoundedChunkWriter {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        if buffer.len() > self.maximum_chunk {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "writer received a buffered document-sized chunk",
            ));
        }
        self.bytes.extend_from_slice(buffer);
        Ok(buffer.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn smiles_writes_molecules_and_repeated_topology_instances_in_authoritative_order() {
    let water = molecule("O");
    let sodium = molecule("[Na+]");
    let mut builder = TopologyBuilder::new();
    let water_definition = builder.add_molecule_definition(&water).unwrap();
    let sodium_definition = builder.add_molecule_definition(&sodium).unwrap();
    builder.add_instance(water_definition).unwrap();
    builder.add_instance(water_definition).unwrap();
    builder.add_instance(sodium_definition).unwrap();
    let topology = builder.build().unwrap();

    assert_eq!(
        smiles::write_molecule(&water, SmilesWriteOptions::default()).unwrap(),
        "O"
    );
    assert_eq!(
        smiles::write_topology(&topology, SmilesWriteOptions::default()).unwrap(),
        "O.O.[Na+]"
    );

    let chiral = molecule("F[C@H](Cl)Br");
    assert_eq!(
        smiles::write_molecule(
            &chiral,
            SmilesWriteOptions {
                mode: SmilesWriteMode::Isomeric,
            },
        )
        .unwrap(),
        smiles::write_isomeric(&chiral).unwrap()
    );
    assert_eq!(
        smiles::write_molecule(
            &molecule("OCC"),
            SmilesWriteOptions {
                mode: SmilesWriteMode::Canonical,
            },
        )
        .unwrap(),
        "CCO"
    );
}

#[test]
fn molfile_model_flattens_disconnected_topology_and_preserves_coordinates() {
    let model = model(
        "CO.O",
        &[[1.25, 2.5, 3.75], [2.0, 2.5, 3.75], [-4.0, 0.5, 8.0]],
    );
    let text = molfile::write_model(&model, MolfileWriteOptions::default()).unwrap();
    assert!(text.contains("V2000"));

    let interpreted = molfile::interpret(&molfile::parse_str(&text).unwrap()).unwrap();
    let round_trip_model = interpreted.into_model();
    assert_eq!(round_trip_model.topology().instance_count(), 2);
    for (index, expected) in [
        Point3::new(1.25, 2.5, 3.75),
        Point3::new(2.0, 2.5, 3.75),
        Point3::new(-4.0, 0.5, 8.0),
    ]
    .into_iter()
    .enumerate()
    {
        let actual = round_trip_model
            .positions()
            .position_at(index)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap();
        assert!((actual.x - expected.x).abs() <= 1e-4);
        assert!((actual.y - expected.y).abs() <= 1e-4);
        assert!((actual.z - expected.z).abs() <= 1e-4);
    }
}

#[test]
fn molfile_version_policy_promotes_counts_and_explicit_v2000_fails() {
    let sodium = molecule("[Na+]");
    let molecules = vec![sodium; 1_000];
    let topology = Topology::from_molecules(&molecules).unwrap();
    let model = Model::new(topology, Positions::zeros(1_000)).unwrap();

    let error = molfile::write_model(
        &model,
        MolfileWriteOptions {
            version: MolfileWriteVersion::V2000,
        },
    )
    .unwrap_err();
    assert_eq!(error.kind(), MolWriteErrorKind::UnsupportedRepresentation);

    let promoted = molfile::write_model(&model, MolfileWriteOptions::default()).unwrap();
    assert!(promoted.contains("V3000"));
    assert!(promoted.contains("M  V30 COUNTS 1000 0 0 0 0"));
    assert_eq!(
        molfile::interpret(&molfile::parse_str(&promoted).unwrap())
            .unwrap()
            .into_model()
            .topology()
            .instance_count(),
        1_000
    );
}

#[test]
fn sdf_writes_independent_models_and_ensemble_members_in_order() {
    let first = model("C", &[[1.0, 0.0, 0.0]]);
    let second = model("O", &[[2.0, 0.0, 0.0]]);
    let text =
        sdf::write_models(&[first.clone(), second.clone()], SdfWriteOptions::default()).unwrap();
    let interpreted = sdf::interpret(&sdf::parse_str(&text).unwrap()).unwrap();
    assert_eq!(interpreted.records().len(), 2);
    assert_eq!(
        interpreted.records()[0]
            .model()
            .topology()
            .atoms()
            .next()
            .unwrap()
            .1
            .element
            .symbol(),
        "C"
    );
    assert_eq!(
        interpreted.records()[1]
            .model()
            .topology()
            .atoms()
            .next()
            .unwrap()
            .1
            .element
            .symbol(),
        "O"
    );

    let shared = first.shared_topology();
    let later = Model::new(
        shared,
        Positions::new(Quantity::new([Point3::new(9.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();
    let ensemble = Ensemble::from_models(&[first, later]).unwrap();
    let ensemble_text = sdf::write_ensemble(&ensemble, SdfWriteOptions::default()).unwrap();
    let records = sdf::interpret(&sdf::parse_str(&ensemble_text).unwrap()).unwrap();
    assert_eq!(records.records().len(), 2);
    assert!(
        (records.records()[0]
            .model()
            .positions()
            .position_at(0)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap()
            .x
            - 1.0)
            .abs()
            <= 1e-9
    );
    assert!(
        (records.records()[1]
            .model()
            .positions()
            .position_at(0)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap()
            .x
            - 9.0)
            .abs()
            <= 1e-9
    );
}

#[test]
fn sdf_explicit_records_preserve_title_and_data_fields() {
    let record = SdfRecordInterpretation::new(
        "named record",
        model("CO.O", &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0]]),
        vec![SdfDataField::new("SOURCE", "round trip")],
    );
    let text = sdf::write_records(&[record], SdfWriteOptions::default()).unwrap();
    let records = sdf::interpret(&sdf::parse_str(&text).unwrap()).unwrap();
    assert_eq!(records.records()[0].title(), "named record");
    assert_eq!(records.records()[0].data_fields()[0].name(), "SOURCE");
    assert_eq!(records.records()[0].data_fields()[0].value(), "round trip");
    assert_eq!(records.records()[0].model().topology().instance_count(), 2);
}

#[test]
fn mmcif_distinguishes_independent_models_from_an_ensemble() {
    let first = model("CO", &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let unrelated = model("N", &[[3.0, 0.0, 0.0]]);
    let independent = mmcif::write_models_with_classifications(
        &[first.clone(), unrelated.clone()],
        &[classifications(&first), classifications(&unrelated)],
        MmcifWriteOptions::default(),
    )
    .unwrap();
    let independent_document = mmcif::parse_str(&independent).unwrap();
    assert_eq!(independent_document.blocks().len(), 2);
    assert!(independent.starts_with("data_model_1\n"));
    assert!(independent.contains("data_model_2\n"));
    assert!(mmcif::interpret_ensemble(
        &independent_document,
        MmcifEnsembleInterpretOptions::default(),
    )
    .is_err());

    let shared = first.shared_topology();
    let second = Model::new(
        Arc::clone(&shared),
        Positions::new(Quantity::new(
            [Point3::new(10.0, 0.0, 0.0), Point3::new(20.0, 0.0, 0.0)],
            ANGSTROM,
        ))
        .unwrap(),
    )
    .unwrap();
    let ensemble = Ensemble::from_models(&[first.clone(), second]).unwrap();
    let ensemble_text = mmcif::write_ensemble_with_classifications(
        &ensemble,
        &classifications(&first),
        MmcifWriteOptions::default(),
    )
    .unwrap();
    let ensemble_document = mmcif::parse_str(&ensemble_text).unwrap();
    assert_eq!(ensemble_document.blocks().len(), 1);
    let parsed =
        mmcif::interpret_ensemble(&ensemble_document, MmcifEnsembleInterpretOptions::default())
            .unwrap();
    assert_eq!(parsed.ensemble().len(), 2);
    assert!(
        (parsed
            .ensemble()
            .member(0)
            .unwrap()
            .positions()
            .position_at(0)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap()
            .x
            - 1.0)
            .abs()
            <= 1e-9
    );
    assert!(
        (parsed
            .ensemble()
            .member(1)
            .unwrap()
            .positions()
            .position_at(0)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap()
            .x
            - 10.0)
            .abs()
            <= 1e-9
    );

    let round_trip =
        mmcif::write_ensemble_interpretation(&parsed, MmcifWriteOptions::default()).unwrap();
    let reparsed = mmcif::interpret_ensemble(
        &mmcif::parse_str(&round_trip).unwrap(),
        MmcifEnsembleInterpretOptions::default(),
    )
    .unwrap();
    assert_eq!(reparsed.ensemble().len(), 2);
}

#[test]
fn mmcif_model_is_one_block_with_coordinate_model_one_and_automatic_classification() {
    let model = model("C", &[[4.0, 5.0, 6.0]]);
    let automatic = mmcif::write_model(&model, MmcifWriteOptions::default()).unwrap();
    assert!(automatic.contains("1 non-polymer"));
    let text = mmcif::write_model_with_classifications(
        &model,
        &classifications(&model),
        MmcifWriteOptions::default(),
    )
    .unwrap();
    assert_eq!(mmcif::parse_str(&text).unwrap().blocks().len(), 1);
    assert!(text.lines().any(|line| line.ends_with(" 1")));
}

#[test]
fn mmcif_export_resolves_reused_atoms_in_hierarchy_order_and_keeps_unassigned_instances() {
    use kekule::core::{Atom, Element, MoleculeEditor};
    use kekule::topology::{AtomSiteMetadata, InstanceAtomId};

    let water = molecule("O");
    let oxygen = water.atom_ids().next().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&water).unwrap();
    let waters = (0..3)
        .map(|_| builder.add_instance(definition).unwrap())
        .collect::<Vec<_>>();
    let mut sodium = Atom::new(Element::from_symbol("Na").unwrap());
    sodium.formal_charge = 1;
    let mut ion = MoleculeEditor::new();
    ion.add_atom(sodium).unwrap();
    builder.add_molecule(&ion.finish().unwrap()).unwrap();
    let hierarchy = builder.hierarchy_mut();
    let chain = hierarchy.add_chain("W", Some("solvent".into())).unwrap();
    // Repeated local AtomIds resolve through their instances; hierarchy order
    // intentionally differs from topology dense order.
    for index in [2, 0, 1] {
        let residue = hierarchy
            .add_residue(chain, "HOH", None, Some((index + 1).to_string()), None)
            .unwrap();
        hierarchy
            .add_atom_site(
                residue,
                InstanceAtomId::new(waters[index], oxygen),
                AtomSiteMetadata {
                    label_atom_id: Some(format!("O{index}")),
                    ..AtomSiteMetadata::default()
                },
            )
            .unwrap();
    }
    let positions = Positions::new(Quantity::new(
        (0..4)
            .map(|i| Point3::new(i as f64 + 1.0, 0.0, 0.0))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))
    .unwrap();
    let mut model = Model::new(builder.build().unwrap(), positions).unwrap();
    for index in 0..3 {
        model
            .set_occupancy(model.atom_ids()[index], Some((index + 1) as f64 / 10.0))
            .unwrap();
    }
    let text = mmcif::write_model(&model, MmcifWriteOptions::default()).unwrap();
    let document = mmcif::parse_str(&text).unwrap();
    let rows = document.blocks()[0].loop_with_tag("_atom_site.id").unwrap();
    assert_eq!(rows.row_count(), 4);
    for (row, index) in [2, 0, 1].into_iter().enumerate() {
        assert_eq!(
            rows.value(row, "_atom_site.label_atom_id").unwrap().text(),
            format!("O{index}")
        );
        assert_eq!(
            rows.value(row, "_atom_site.auth_seq_id").unwrap().text(),
            (index + 1).to_string()
        );
        assert_eq!(
            rows.value(row, "_atom_site.Cartn_x")
                .unwrap()
                .text()
                .parse::<f64>()
                .unwrap(),
            index as f64 + 1.0
        );
        assert_eq!(
            rows.value(row, "_atom_site.occupancy")
                .unwrap()
                .text()
                .parse::<f64>()
                .unwrap(),
            (index + 1) as f64 / 10.0
        );
    }
    assert_eq!(
        rows.value(3, "_atom_site.type_symbol").unwrap().text(),
        "Na"
    );
    assert_eq!(
        rows.value(3, "_atom_site.Cartn_x")
            .unwrap()
            .text()
            .parse::<f64>()
            .unwrap(),
        4.0
    );
}

#[test]
fn mmcif_export_still_rejects_a_partially_assigned_molecule() {
    use kekule::topology::{AtomSiteMetadata, InstanceAtomId};

    let molecule = molecule("CO");
    let atoms = molecule.atom_ids().collect::<Vec<_>>();
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&molecule).unwrap();
    let hierarchy = builder.hierarchy_mut();
    let chain = hierarchy.add_chain("A", None).unwrap();
    let residue = hierarchy
        .add_residue(chain, "UNL", None, None, None)
        .unwrap();
    // Put only the second atom in hierarchy so testing the first atom alone
    // cannot incorrectly choose the hierarchy-free output path.
    hierarchy
        .add_atom_site(
            residue,
            InstanceAtomId::new(instance, atoms[1]),
            AtomSiteMetadata::default(),
        )
        .unwrap();
    let model = Model::new(builder.build().unwrap(), Positions::zeros(2)).unwrap();
    assert_eq!(
        mmcif::write_model(&model, MmcifWriteOptions::default()),
        Err(mmcif::MmcifWriteError::MissingAtomSite(
            InstanceAtomId::new(instance, atoms[0])
        ))
    );
}

#[test]
fn mmcif_streaming_matches_string_output_for_models_and_ensemble() {
    let first = model("CO", &[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let unrelated = model("N", &[[3.0, 0.0, 0.0]]);
    let models = [first.clone(), unrelated.clone()];
    let model_classifications = [classifications(&first), classifications(&unrelated)];
    let expected_models = mmcif::write_models_with_classifications(
        &models,
        &model_classifications,
        MmcifWriteOptions::default(),
    )
    .unwrap();
    let mut streamed_models = BoundedChunkWriter {
        bytes: Vec::new(),
        maximum_chunk: 64,
    };
    mmcif::write_models_with_classifications_to(
        &mut streamed_models,
        &models,
        &model_classifications,
        MmcifWriteOptions::default(),
    )
    .unwrap();
    assert_eq!(streamed_models.bytes, expected_models.as_bytes());
    assert_eq!(
        mmcif::parse_str(std::str::from_utf8(&streamed_models.bytes).unwrap())
            .unwrap()
            .blocks()
            .len(),
        2
    );

    let second = Model::new(
        first.shared_topology(),
        Positions::new(Quantity::new(
            [Point3::new(10.0, 0.0, 0.0), Point3::new(20.0, 0.0, 0.0)],
            ANGSTROM,
        ))
        .unwrap(),
    )
    .unwrap();
    let ensemble = Ensemble::from_models(&[first.clone(), second]).unwrap();
    let expected_ensemble = mmcif::write_ensemble_with_classifications(
        &ensemble,
        &classifications(&first),
        MmcifWriteOptions::default(),
    )
    .unwrap();
    let mut streamed_ensemble = BoundedChunkWriter {
        bytes: Vec::new(),
        maximum_chunk: 64,
    };
    mmcif::write_ensemble_with_classifications_to(
        &mut streamed_ensemble,
        &ensemble,
        &classifications(&first),
        MmcifWriteOptions::default(),
    )
    .unwrap();
    assert_eq!(streamed_ensemble.bytes, expected_ensemble.as_bytes());

    let document =
        mmcif::parse_str(std::str::from_utf8(&streamed_ensemble.bytes).unwrap()).unwrap();
    assert_eq!(document.blocks().len(), 1);
    assert_eq!(
        mmcif::interpret_ensemble(&document, MmcifEnsembleInterpretOptions::default())
            .unwrap()
            .ensemble()
            .len(),
        2
    );
}

#[test]
fn mmcif_report_count_mismatch_is_context_neutral() {
    let first = model("C", &[[1.0, 0.0, 0.0]]);
    let second = Model::new(
        first.shared_topology(),
        Positions::new(Quantity::new([Point3::new(2.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();

    let models_error = mmcif::write_models_with_reports(
        &[first.clone(), second.clone()],
        &[],
        MmcifWriteOptions::default(),
    )
    .unwrap_err();
    assert!(matches!(
        &models_error,
        mmcif::MmcifWriteError::ReportCountMismatch {
            expected: 2,
            actual: 0
        }
    ));
    assert_eq!(
        models_error.to_string(),
        "mmCIF writing requires 2 interpretation reports, but received 0"
    );

    let ensemble = Ensemble::from_models(&[first, second]).unwrap();
    let ensemble_error =
        mmcif::write_ensemble_with_reports(&ensemble, &[], MmcifWriteOptions::default())
            .unwrap_err();
    assert!(matches!(
        &ensemble_error,
        mmcif::MmcifWriteError::ReportCountMismatch {
            expected: 2,
            actual: 0
        }
    ));
    assert_eq!(ensemble_error.to_string(), models_error.to_string());
}

struct FailingWriter;

impl Write for FailingWriter {
    fn write(&mut self, _buffer: &[u8]) -> io::Result<usize> {
        Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed"))
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn streaming_writers_report_structured_io_failures() {
    let model = model("C", &[[0.0, 0.0, 0.0]]);
    assert!(matches!(
        molfile::write_model_to(
            &mut FailingWriter,
            &model,
            MolfileWriteOptions::default()
        ),
        Err(error) if error.kind() == MolWriteErrorKind::Io(io::ErrorKind::BrokenPipe)
    ));
    assert!(matches!(
        sdf::write_model_to(&mut FailingWriter, &model, SdfWriteOptions::default()),
        Err(sdf::SdfWriteError::Io {
            kind: io::ErrorKind::BrokenPipe,
            ..
        })
    ));
    assert!(matches!(
        mmcif::write_model_with_classifications_to(
            &mut FailingWriter,
            &model,
            &classifications(&model),
            MmcifWriteOptions::default(),
        ),
        Err(mmcif::MmcifWriteError::Io {
            kind: io::ErrorKind::BrokenPipe,
            ..
        })
    ));
}
