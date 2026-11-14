use super::*;

/// A deterministic linear congruential generator.
struct Lcg(u64);

impl Lcg {
    fn below(&mut self, bound: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) % bound
    }
}

type Contacts = (Vec<XyzInferredBond>, Vec<XyzInterpretIssue>);

/// Every atom pair measured: each atom meets the later atoms within
/// bonding distance in atom-line order.
fn brute(
    elements: &[Element],
    points: &[Point3],
    max_contacts: usize,
) -> Result<Contacts, (LimitExceeded, usize)> {
    let radius = |atom: usize| covalent_bonding_radius(elements[atom]).0;
    let mut contacts = vec![0; points.len()];
    let (mut bonds, mut close) = (Vec::new(), Vec::new());
    for left in 0..points.len() {
        for right in left + 1..points.len() {
            let (point, other) = (points[left], points[right]);
            let distance_squared = (point.x - other.x).powi(2)
                + (point.y - other.y).powi(2)
                + (point.z - other.z).powi(2);
            let cutoff = radius(left) + radius(right) + COVALENT_BOND_TOLERANCE_ANGSTROM;
            if distance_squared > cutoff * cutoff {
                continue;
            }
            for atom in [left, right] {
                contacts[atom] += 1;
                if contacts[atom] > max_contacts {
                    return Err((
                        LimitExceeded::new(Resource::CandidatePairs, max_contacts, contacts[atom]),
                        atom,
                    ));
                }
            }
            let atoms = [dense_index(left), dense_index(right)];
            let distance = distance_squared.sqrt();
            if distance_squared <= MIN_COVALENT_BOND_DISTANCE_ANGSTROM.powi(2) {
                close.push(XyzInterpretIssue::CloseContact {
                    atoms,
                    distance: Quantity::new(distance, ANGSTROM),
                });
            } else {
                let order = if is_metal(elements[left]) != is_metal(elements[right]) {
                    BondOrder::Dative
                } else {
                    BondOrder::Single
                };
                bonds.push(XyzInferredBond {
                    atoms,
                    order,
                    distance,
                    cutoff,
                });
            }
        }
    }
    Ok((bonds, close))
}

/// The grid finds the contacts every pair gives, and with a low limit
/// fails at the same atom; atoms sit on cell boundaries and exactly at
/// bonding distance.
#[test]
fn contacts_match_every_pair() {
    let palette = ["H", "C", "O", "Fe", "Cl"]
        .map(|symbol| Element::from_symbol(symbol).expect("element symbol"));
    for seed in 0..12 {
        let mut random = Lcg(seed);
        let count = 300;
        let spread = [6, 12, 40][seed as usize % 3];
        let elements = (0..count)
            .map(|_| palette[random.below(palette.len() as u64) as usize])
            .collect::<Vec<_>>();
        let mut points = Vec::<Point3>::new();
        for index in 0..count {
            let point = if index % 7 == 3 {
                // Exactly at bonding distance from the previous atom.
                let previous = points[index - 1];
                let cutoff = covalent_bonding_radius(elements[index - 1]).0
                    + covalent_bonding_radius(elements[index]).0
                    + COVALENT_BOND_TOLERANCE_ANGSTROM;
                Point3::new(previous.x, previous.y + cutoff, previous.z)
            } else if index % 7 == 5 {
                // On half-ångström steps, often on cell boundaries.
                let mut axis = || random.below(spread * 2) as f64 * 0.5 - spread as f64 * 0.5;
                Point3::new(axis(), axis(), axis())
            } else {
                let mut axis = || random.below(spread * 1000) as f64 * 1e-3 - spread as f64 * 0.5;
                Point3::new(axis(), axis(), axis())
            };
            points.push(point);
        }
        let (mut bonds, mut close) = brute(&elements, &points, usize::MAX).unwrap();
        assert!(!bonds.is_empty(), "seed {seed}");
        bonds.sort_by_key(|bond| bond.atoms);
        close.sort_by_key(|issue| match issue {
            XyzInterpretIssue::CloseContact { atoms, .. } => *atoms,
            _ => unreachable!("only close contacts"),
        });
        assert_eq!(
            infer_contacts(&elements, &points, usize::MAX).unwrap(),
            (bonds, close),
            "seed {seed}"
        );
        for limit in 0..6 {
            let expected = brute(&elements, &points, limit).map(|_| ());
            let found = infer_contacts(&elements, &points, limit).map(|_| ());
            assert_eq!(found, expected, "seed {seed}, limit {limit}");
        }
    }
}
