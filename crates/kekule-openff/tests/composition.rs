use kekule::{
    smiles,
    topology::{Topology, TopologyBuilder},
    units::{ELEMENTARY_CHARGE, NANOMETER},
};
use kekule_openff::{ChargeSource, ElectrostaticsMethod, ForceField, VdwMethod};
use std::sync::Arc;

// Deliberately simple numerical parameters isolate assignment contracts.
const WATER: &str = r#"<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
<Bonds version="0.4"><Bond smirks="[#8:1]-[#1:2]" id="water-bond" length="0.1*nanometer" k="100*kilojoule_per_mole/nanometer**2"/></Bonds>
<Angles version="0.3"><Angle smirks="[#1:1]-[#8:2]-[#1:3]" id="water-angle" angle="90*degree" k="100*kilojoule_per_mole/radian**2"/></Angles>
<ProperTorsions version="0.4"/>
<Constraints version="0.3"><Constraint smirks="[#1:1]-[#8]-[#1:2]" id="water-hh" distance="0.15*nanometer"/></Constraints>
<vdW version="0.4"><Atom smirks="[*:1]" id="water-vdw" sigma="0.3*nanometer" epsilon="0.1*kilojoule_per_mole"/></vdW>
<Electrostatics version="0.4" scale14="0.8333333333"/>
<LibraryCharges version="0.3"><LibraryCharge smirks="[#8:1](-[#1:2])-[#1:3]" id="water-charge" charge1="-0.8*elementary_charge" charge2="0.4*elementary_charge" charge3="0.4*elementary_charge"/></LibraryCharges>
</SMIRNOFF>"#;

fn molecule(text: &str) -> kekule::core::Molecule {
    let mut m = smiles::to_molecules(text).unwrap().remove(0);
    m.perceive().unwrap();
    m.add_hydrogens().unwrap();
    m
}

#[test]
fn library_only_repeated_definitions_retain_exact_binding_and_nonbonded_constraint() {
    let ff = ForceField::from_offxml(WATER).unwrap();
    let m = molecule("O");
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(m.clone()).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    let topology = Arc::new(builder.build().unwrap());
    let p = ff.parameterize_without_nagl(Arc::clone(&topology)).unwrap();
    assert!(Arc::ptr_eq(p.topology(), &topology));
    assert_eq!(
        p.charges().value_in(ELEMENTARY_CHARGE).unwrap(),
        [-0.8, 0.4, 0.4, -0.8, 0.4, 0.4]
    );
    assert_eq!(p.constraints().len(), 2);
    assert_ne!(p.constraints()[0].atoms, p.constraints()[1].atoms);
    for c in p.constraints() {
        assert_eq!(c.parameter.distance.value_in(NANOMETER).unwrap(), 0.15);
    }
    assert!(p.charge_sources().iter().all(|s| matches!(s, ChargeSource::Library { parameter_ids } if parameter_ids == &["water-charge"])));
}

#[test]
fn composition_is_ordered_and_keeps_rosemary_model_and_unrelated_parameters() {
    let mut ff = ForceField::rosemary().unwrap();
    let model = ff.charge_model().cloned();
    let organic = molecule("CC");
    let original = ff.label_molecule(&organic).unwrap();
    // Restrict the synthetic vdW rule to water so unrelated atoms are untouched.
    let xml = WATER.replace("[*:1]", "[$([#8X2H2]),$([#1]-[#8X2H2]):1]");
    let water = ForceField::from_offxml(&xml).unwrap();
    ff.append(&water).unwrap();
    assert_eq!(ff.charge_model(), model.as_ref());
    assert_eq!(ff.label_molecule(&organic).unwrap().bonds, original.bonds);
    assert_eq!(ff.label_molecule(&organic).unwrap().vdw, original.vdw);
    let topology = Arc::new(Topology::from_molecule((molecule("O")).clone()).unwrap());
    let p = ff.parameterize_without_nagl(Arc::clone(&topology)).unwrap();
    assert!(p
        .bonds()
        .iter()
        .all(|t| t.parameter.source.id == "water-bond"));
    assert_eq!(p.charges().value()[0], -0.8);
    let newer = ForceField::from_offxml(
        &xml.replace("water-charge", "new-charge")
            .replace("-0.8*", "-0.6*")
            .replace("0.4*", "0.3*"),
    )
    .unwrap();
    ff.append(&newer).unwrap();
    assert_eq!(
        ff.parameterize_without_nagl(topology)
            .unwrap()
            .charges()
            .value()[0],
        -0.6
    );
    let mut reversed = water;
    reversed.append(&ForceField::rosemary().unwrap()).unwrap();
    assert_eq!(reversed.charge_model(), model.as_ref());
    assert!(reversed
        .label_molecule(&molecule("O"))
        .unwrap()
        .bonds
        .values()
        .all(|p| p.id != "water-bond"));
}

#[test]
fn failed_composition_is_atomic_for_settings_and_model_conflicts() {
    let mut ff = ForceField::rosemary().unwrap();
    let original = ff.label_molecule(&molecule("O")).unwrap();
    let incompatible = ForceField::from_offxml(
        &WATER.replace("<vdW version", "<vdW cutoff=\"1.2*nanometer\" version"),
    )
    .unwrap();
    assert!(ff
        .append(&incompatible)
        .unwrap_err()
        .to_string()
        .contains("nonbonded settings"));
    let rosemary = include_str!("../data/rosemary.offxml");
    let incompatible =
        ForceField::from_offxml(&rosemary.replace("openff-gnn-am1bcc-1.0.0.pt", "different.pt"))
            .unwrap();
    assert!(ff
        .append(&incompatible)
        .unwrap_err()
        .to_string()
        .contains("NAGL models"));
    assert_eq!(
        ff.label_molecule(&molecule("O")).unwrap().bonds,
        original.bonds
    );
    assert_eq!(
        ff.charge_model().unwrap().model_file(),
        "openff-gnn-am1bcc-1.0.0.pt"
    );
}

#[test]
fn composition_accepts_length_unit_roundoff_but_rejects_changed_scales() {
    let mut ff = ForceField::rosemary().unwrap();
    let xml = WATER
        .replace(
            "<vdW version",
            "<vdW cutoff=\"0.9*nanometer\" switch_width=\"0.1*nanometer\" version",
        )
        .replace(
            "<Electrostatics version",
            "<Electrostatics cutoff=\"0.9*nanometer\" version",
        );
    let settings = ff.nonbonded_settings().clone();
    ff.append(&ForceField::from_offxml(&xml).unwrap()).unwrap();
    assert_eq!(ff.nonbonded_settings(), &settings);
    assert!(ff
        .append(
            &ForceField::from_offxml(&xml.replace("0.8333333333", "0.8333333333333334")).unwrap()
        )
        .is_err());
}

#[test]
fn incomplete_or_nonconserving_charges_and_implicit_nonbonded_distances_fail() {
    let topology = Arc::new(Topology::from_molecule((molecule("O")).clone()).unwrap());
    for xml in [
        WATER.replace(
            "charge1=\"-0.8*elementary_charge\"",
            "charge1=\"-0.7*elementary_charge\"",
        ),
        WATER.replace("[#8:1](-[#1:2])-[#1:3]", "[#6:1](-[#1:2])-[#1:3]"),
        WATER.replace(" distance=\"0.15*nanometer\"", ""),
    ] {
        let ff = ForceField::from_offxml(&xml).unwrap();
        assert!(ff.parameterize_without_nagl(Arc::clone(&topology)).is_err());
    }
    let ff = ForceField::rosemary().unwrap();
    let organic = Arc::new(Topology::from_molecule((molecule("CC")).clone()).unwrap());
    assert!(ff
        .parameterize_without_nagl(organic)
        .unwrap_err()
        .to_string()
        .contains("incomplete library charges"));
}

fn two_waters() -> Arc<Topology> {
    let m = molecule("O");
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(m.clone()).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    Arc::new(builder.build().unwrap())
}

#[test]
fn interactions_from_one_rule_share_one_parameter_allocation() {
    let p = ForceField::from_offxml(WATER)
        .unwrap()
        .parameterize_without_nagl(two_waters())
        .unwrap();
    let first = &p.bonds()[0];
    assert!(p
        .bonds()
        .iter()
        .all(|b| Arc::ptr_eq(&b.parameter, &first.parameter)));
    assert!(p
        .bonds()
        .iter()
        .any(|b| b.atoms[0].molecule() != first.atoms[0].molecule()));
    assert!(p
        .angles()
        .iter()
        .all(|a| Arc::ptr_eq(&a.parameter, &p.angles()[0].parameter)));
    assert!(p
        .constraints()
        .iter()
        .all(|c| Arc::ptr_eq(&c.parameter, &p.constraints()[0].parameter)));
    assert!(p.vdw().iter().all(|v| Arc::ptr_eq(v, &p.vdw()[0])));
}

#[test]
fn nonbonded_methods_are_typed_and_keep_smirnoff_spelling() {
    let settings = |xml: &str| {
        ForceField::from_offxml(xml)
            .unwrap()
            .parameterize_without_nagl(two_waters())
            .unwrap()
            .nonbonded_settings()
            .clone()
    };
    let defaults = settings(WATER);
    assert_eq!(defaults.vdw_periodic_method, VdwMethod::Cutoff);
    assert_eq!(defaults.vdw_nonperiodic_method, VdwMethod::NoCutoff);
    assert_eq!(
        defaults.electrostatics_periodic_method,
        ElectrostaticsMethod::Ewald3DConductingBoundary
    );
    assert_eq!(
        defaults.electrostatics_nonperiodic_method,
        ElectrostaticsMethod::Coulomb
    );
    assert_eq!(VdwMethod::NoCutoff.smirnoff_name(), "no-cutoff");
    assert_eq!(
        ElectrostaticsMethod::Ewald3DConductingBoundary.smirnoff_name(),
        "Ewald3D-ConductingBoundary"
    );
    let changed = settings(
        &WATER
            .replace(
                r#"<vdW version="0.4">"#,
                r#"<vdW version="0.4" periodic_method="no-cutoff" nonperiodic_method="cutoff">"#,
            )
            .replace(
                r#"<Electrostatics version="0.4""#,
                r#"<Electrostatics version="0.4" periodic_potential="Coulomb""#,
            ),
    );
    assert_eq!(changed.vdw_periodic_method, VdwMethod::NoCutoff);
    assert_eq!(changed.vdw_nonperiodic_method, VdwMethod::Cutoff);
    assert_eq!(
        changed.electrostatics_periodic_method,
        ElectrostaticsMethod::Coulomb
    );
    let pme = settings(&WATER.replace(
        r#"<Electrostatics version="0.4""#,
        r#"<Electrostatics version="0.4" periodic_potential="PME""#,
    ));
    assert_eq!(
        pme.electrostatics_periodic_method,
        ElectrostaticsMethod::Ewald3DConductingBoundary
    );
}
