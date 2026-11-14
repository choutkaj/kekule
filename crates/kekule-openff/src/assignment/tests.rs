use super::*;
use kekule::diagnostics::{LimitExceeded, Resource};
use kekule::substructure::SubstructureMatchWork;
use kekule::{smiles, topology::TopologyBuilder};

#[test]
fn exceeded_match_limits_are_resource_limits_with_the_wrapped_cause() {
    let limit = LimitExceeded::new(Resource::SearchStates, 10, 11);
    let error = matching(
        "parameter \"b1\"",
        TaggedMatchError::Match(SubstructureMatchError::LimitExceeded {
            limit,
            work: SubstructureMatchWork::default(),
        }),
    );
    assert_eq!(error.kind(), ErrorKind::ResourceLimit);
    assert_eq!(error.limit_exceeded(), Some(limit));
    assert_eq!(
        error.to_string(),
        "Cannot match the SMIRKS pattern of parameter \"b1\"."
    );
    assert!(matches!(
        std::error::Error::source(&error)
            .unwrap()
            .downcast_ref::<TaggedMatchError>(),
        Some(TaggedMatchError::Match(
            SubstructureMatchError::LimitExceeded { .. }
        ))
    ));

    let error = matching("parameter \"b1\"", TaggedMatchError::ZeroTag);
    assert_eq!(error.kind(), ErrorKind::ForceField);
    assert_eq!(error.limit_exceeded(), None);
    assert!(matches!(
        std::error::Error::source(&error)
            .unwrap()
            .downcast_ref::<TaggedMatchError>(),
        Some(TaggedMatchError::ZeroTag)
    ));
}

#[test]
fn fourteen_scaling_uses_shortest_bond_distance() {
    let mut molecule = smiles::to_molecules("CCCC").unwrap().remove(0);
    molecule.perceive().unwrap();
    molecule.add_hydrogens(Default::default()).unwrap();
    let carbons = molecule
        .atoms()
        .filter(|(_, a)| a.element().atomic_number() == 6)
        .map(|(a, _)| a)
        .collect::<Vec<_>>();
    let topology = Arc::new(Topology::from_molecule(molecule.clone()).unwrap());
    let ff = ForceField::rosemary().unwrap();
    let p = ff
        .parameterize_with(topology, |m| {
            Ok(ChargeAssignment {
                charges: Quantity::new(vec![0.0; m.atom_count()], ELEMENTARY_CHARGE),
                source: ChargeSource::Library {
                    parameter_ids: vec![],
                },
            })
        })
        .unwrap();
    let pair = p
        .pair_exceptions()
        .iter()
        .find(|p| p.atoms.map(|a| a.atom()) == [carbons[0], carbons[3]])
        .unwrap();
    assert_eq!(pair.vdw_scale, 0.5);
    assert_eq!(pair.electrostatics_scale, 0.8333333333);
}

#[test]
fn reused_definitions_snapshot_trefoils_constraints_and_exceptions() {
    let molecule = smiles::to_molecules("[H]C([H])=O").unwrap().remove(0);
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(molecule.clone()).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    let topology = Arc::new(builder.build().unwrap());
    let ff = ForceField::rosemary().unwrap();
    let mut calls = 0;
    let p = ff
        .parameterize_with(topology.clone(), |m| {
            calls += 1;
            Ok(ChargeAssignment {
                charges: Quantity::new(vec![0.0; m.atom_count()], ELEMENTARY_CHARGE),
                source: ChargeSource::Library {
                    parameter_ids: vec![],
                },
            })
        })
        .unwrap();
    assert_eq!(calls, 1);
    assert!(Arc::ptr_eq(&p.shared_topology(), &topology));
    assert_eq!(p.charges().value(), &[0.0; 8]);
    assert_eq!(p.bonds().len(), 6);
    assert_eq!(p.angles().len(), 6);
    assert_eq!(p.improper_torsions().len(), 6);
    for t in p.improper_torsions() {
        let center = topology.atom(t.atoms[0]).unwrap();
        assert_eq!(center.element().atomic_number(), 6);
        assert!(t.parameter.terms.iter().all(|term| term.idivf == 3.0));
    }
    assert_eq!(p.constraints().len(), 4);
    for c in p.constraints() {
        let b = p.bonds().iter().find(|b| b.atoms == c.atoms).unwrap();
        assert_eq!(c.parameter.distance, b.parameter.length);
    }
    // All six pairs within each four-atom star are at distance one or two;
    // instance boundaries must never introduce exceptions.
    assert_eq!(p.pair_exceptions().len(), 12);
    assert!(p
        .pair_exceptions()
        .iter()
        .all(|e| e.vdw_scale == 0.0 && e.electrostatics_scale == 0.0));
}

#[test]
fn complete_library_charges_precede_inference_partial_coverage_falls_back() {
    let ff = ForceField::rosemary().unwrap();
    let ion = crate::explicit(&smiles::to_molecules("[Na+]").unwrap().remove(0)).unwrap();
    let charges = ff
        .charges(&ion, |_| panic!("library ion must not invoke NAGL"))
        .unwrap();
    assert_eq!(charges.charges.value(), &[1.0]);
    assert!(matches!(charges.source, ChargeSource::Library { .. }));
    let methanol = crate::explicit(
        &smiles::to_molecules("[H]C([H])([H])O[H]")
            .unwrap()
            .remove(0),
    )
    .unwrap();
    let mut modified = ff;
    let h = crate::offxml::LibraryCharge {
        source: ParameterIdentity {
            id: "partial".into(),
            smirks: "[#1:1]".into(),
        },
        charges: vec![0.0],
    };
    modified.library.push(Rule {
        query: kekule::smarts::parse_str(&h.source.smirks, Default::default()).unwrap(),
        parameter: h,
    });
    assert_eq!(
        modified
            .charges(&methanol, |_| Err(Error::new(
                ErrorKind::Charges,
                "fallback reached"
            )))
            .unwrap_err()
            .to_string(),
        "fallback reached"
    );
}
