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
