use super::*;
use crate::algorithms::allowed_valences_with_charge;

#[test]
fn every_builder_state_is_a_valence_the_model_accepts() {
    for (number, charge, valences) in VALENCES {
        let charges = CHARGES
            .iter()
            .find(|(element, _)| element == number)
            .map(|(_, charges)| *charges)
            .unwrap();
        assert!(charges.contains(charge), "{number} {charge}");
        let accepted = allowed_valences_with_charge(*number, *charge).unwrap();
        for valence in *valences {
            assert!(accepted.contains(valence), "{number} {charge} {valence}");
        }
    }
    for (number, charges) in CHARGES {
        assert!(charges.windows(2).all(|pair| pair[0] < pair[1]));
        for charge in *charges {
            assert!(VALENCES
                .iter()
                .any(|(element, state, _)| element == number && state == charge));
        }
    }
}

#[test]
fn errors_name_the_state_and_what_it_takes() {
    let oxygen = builder_hydrogen_count(Element::O, 2, 0).unwrap_err();
    assert_eq!(
        oxygen.to_string(),
        "The atom is O with charge +2. O takes charge -1, 0 or +1."
    );
    let carbon = builder_hydrogen_count(Element::C, 0, 5).unwrap_err();
    assert_eq!(
        carbon.to_string(),
        "The atom has valence 5. C takes valence 4."
    );
    let sulfur = builder_hydrogen_count(Element::S, 1, 4).unwrap_err();
    assert_eq!(
        sulfur.to_string(),
        "The atom has valence 4. S with charge +1 takes valence 3."
    );
    let sulfur = builder_hydrogen_count(Element::S, 0, 7).unwrap_err();
    assert_eq!(
        sulfur.to_string(),
        "The atom has valence 7. S takes valence 2, 4 or 6."
    );
    let iron = builder_hydrogen_count(Element::from_symbol("Fe").unwrap(), 0, 0).unwrap_err();
    assert_eq!(
        iron.to_string(),
        "The atom is Fe. Fe has no builder valences."
    );
}
