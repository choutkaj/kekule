use super::*;
use crate::structure::{Model, Positions};
use crate::topology::{AtomSelection, InstanceAtomId, Topology, TopologyAtomIndex};
use crate::units::ANGSTROM;
use crate::{mmcif, pdb, smiles};

/// One PDB coordinate record.
fn atom(serial: usize, name: &str, residue: &str, chain: char, number: i32, code: char) -> String {
    let field = if name.len() < 4 {
        format!(" {name:<3}")
    } else {
        name.to_owned()
    };
    let x = serial as f64 * 1.3;
    let y = (serial * serial % 7) as f64 * 0.9;
    let z = (serial % 3) as f64 * 1.1;
    format!(
        "ATOM  {serial:>5} {field} {residue} {chain}{number:>4}{code}   {x:8.3}{y:8.3}{z:8.3}  1.00  0.00           {}\n",
        &name[..1],
    )
}

fn read_pdb(records: &[(&str, &str, char, i32, char)]) -> Model {
    let mut text: String = records
        .iter()
        .enumerate()
        .map(|(serial, &(name, residue, chain, number, code))| {
            atom(serial + 1, name, residue, chain, number, code)
        })
        .collect();
    text.push_str("END\n");
    pdb::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .into_model()
}

fn all(model: &Model) -> AtomSelection {
    AtomSelection::all(model.shared_topology())
}

fn ids(topology: &Topology, indices: &[usize]) -> Vec<InstanceAtomId> {
    indices
        .iter()
        .map(|&index| topology.atom_ids()[index])
        .collect()
}

fn paired_indices(paired: &LabelCorrespondence) -> Vec<(usize, usize)> {
    paired
        .correspondence()
        .index_pairs()
        .iter()
        .map(|(moving, reference)| (moving.index(), reference.index()))
        .collect()
}

fn unpaired_indices(paired: &LabelCorrespondence, side: CorrespondenceSide) -> Vec<usize> {
    paired
        .unpaired(side)
        .indices()
        .iter()
        .map(|index| index.index())
        .collect()
}

#[test]
fn author_labels_pair_atoms_of_different_structures_and_report_the_rest() {
    // The moving structure lists the residues in another order, mutates
    // residue 2, and lacks the reference's chain B.
    let moving = read_pdb(&[
        ("N", "ALA", 'A', 2, ' '),
        ("CA", "ALA", 'A', 2, ' '),
        ("CB", "ALA", 'A', 2, ' '),
        ("N", "GLY", 'A', 1, ' '),
        ("CA", "GLY", 'A', 1, ' '),
        ("OXT", "GLY", 'A', 1, ' '),
    ]);
    let reference = read_pdb(&[
        ("N", "GLY", 'A', 1, ' '),
        ("CA", "GLY", 'A', 1, ' '),
        ("N", "SER", 'A', 2, ' '),
        ("CA", "SER", 'A', 2, ' '),
        ("CB", "SER", 'A', 2, ' '),
        ("OG", "SER", 'A', 2, ' '),
        ("N", "GLY", 'B', 1, ' '),
    ]);
    let paired =
        AtomCorrespondence::from_labels(&all(&moving), &all(&reference), LabelScheme::Author)
            .unwrap();
    // Pairs follow moving dense order; residue names are not compared.
    assert_eq!(
        paired_indices(&paired),
        [(0, 2), (1, 3), (2, 4), (3, 0), (4, 1)]
    );
    assert_eq!(unpaired_indices(&paired, CorrespondenceSide::Moving), [5]);
    assert_eq!(
        unpaired_indices(&paired, CorrespondenceSide::Reference),
        [5, 6]
    );
    assert!(paired
        .unpaired_moving()
        .topology()
        .shares_layout(moving.topology()));
    assert!(paired
        .unpaired_reference()
        .topology()
        .shares_layout(reference.topology()));
}

#[test]
fn a_label_correspondence_feeds_the_fits() {
    let moving = read_pdb(&[
        ("N", "GLY", 'A', 1, ' '),
        ("CA", "GLY", 'A', 1, ' '),
        ("C", "GLY", 'A', 1, ' '),
        ("O", "GLY", 'A', 1, ' '),
    ]);
    // The reference holds the same atoms in reverse order, turned and moved.
    let records = [
        ("O", "GLY", 'A', 1, ' '),
        ("C", "GLY", 'A', 1, ' '),
        ("CA", "GLY", 'A', 1, ' '),
        ("N", "GLY", 'A', 1, ' '),
    ];
    let reference = read_pdb(&records);
    let source = moving.positions().values().into_value();
    let turned: Vec<_> = source
        .iter()
        .rev()
        .map(|p| Point3::new(-p.y + 2., p.x - 3., p.z + 1.))
        .collect();
    let reference = Model::new(
        reference.shared_topology(),
        Positions::new(Quantity::new(turned, CANONICAL_LENGTH_UNIT)).unwrap(),
    )
    .unwrap();
    let paired =
        AtomCorrespondence::from_labels(&all(&moving), &all(&reference), LabelScheme::Author)
            .unwrap();
    assert_eq!(paired_indices(&paired), [(0, 3), (1, 2), (2, 1), (3, 0)]);
    let result = fit(
        moving.as_model_view(),
        reference.as_model_view(),
        &paired,
        AlignmentOptions::default(),
    )
    .unwrap();
    assert!(result.rmsd().value_in(ANGSTROM).unwrap() < 1e-9);
    let correspondence: AtomCorrespondence = paired.into();
    assert_eq!(correspondence.len(), 4);
}

#[test]
fn insertion_codes_tell_residues_of_one_number_apart() {
    let moving = read_pdb(&[("CA", "GLY", 'A', 10, 'A'), ("CA", "GLY", 'A', 10, ' ')]);
    let reference = read_pdb(&[
        ("CA", "GLY", 'A', 10, ' '),
        ("CA", "GLY", 'A', 10, 'A'),
        ("CA", "GLY", 'A', 10, 'B'),
    ]);
    let paired =
        AtomCorrespondence::from_labels(&all(&moving), &all(&reference), LabelScheme::Author)
            .unwrap();
    assert_eq!(paired_indices(&paired), [(0, 1), (1, 0)]);
    assert_eq!(
        unpaired_indices(&paired, CorrespondenceSide::Reference),
        [2]
    );
}

#[test]
fn selections_choose_the_atoms_compared() {
    let moving = read_pdb(&[
        ("N", "GLY", 'A', 1, ' '),
        ("CA", "GLY", 'A', 1, ' '),
        ("C", "GLY", 'A', 1, ' '),
    ]);
    let reference = moving.clone();
    let moving_atoms =
        AtomSelection::from_indices(moving.shared_topology(), [1, 2].map(TopologyAtomIndex::new))
            .unwrap();
    let reference_atoms = AtomSelection::from_indices(
        reference.shared_topology(),
        [0, 1].map(TopologyAtomIndex::new),
    )
    .unwrap();
    let paired =
        AtomCorrespondence::from_labels(&moving_atoms, &reference_atoms, LabelScheme::Author)
            .unwrap();
    assert_eq!(paired_indices(&paired), [(1, 1)]);
    // Only selected atoms are reported unpaired.
    assert_eq!(unpaired_indices(&paired, CorrespondenceSide::Moving), [2]);
    assert_eq!(
        unpaired_indices(&paired, CorrespondenceSide::Reference),
        [0]
    );
}

#[test]
fn labels_naming_several_atoms_are_refused_by_name() {
    // A residue number recurring after another residue starts a new residue
    // with the same labels.
    let moving = read_pdb(&[
        ("CA", "GLY", 'A', 1, ' '),
        ("N", "GLY", 'A', 1, ' '),
        ("CA", "GLY", 'A', 2, ' '),
        ("CA", "GLY", 'A', 1, ' '),
        ("N", "GLY", 'A', 1, ' '),
    ]);
    let reference = read_pdb(&[("CA", "GLY", 'A', 1, ' ')]);
    let error =
        AtomCorrespondence::from_labels(&all(&moving), &all(&reference), LabelScheme::Author)
            .unwrap_err();
    let AtomCorrespondenceError::AmbiguousLabels { side, labels } = &error else {
        panic!("{error:?}");
    };
    assert_eq!(*side, CorrespondenceSide::Moving);
    assert_eq!(labels.len(), 2);
    let topology = moving.topology();
    assert_eq!(labels[0].label().chain(), "A");
    assert_eq!(labels[0].label().residue_number(), 1);
    assert_eq!(labels[0].label().insertion_code(), None);
    assert_eq!(labels[0].label().atom_name(), "CA");
    assert_eq!(labels[0].atoms(), ids(topology, &[0, 3]));
    assert_eq!(labels[1].label().atom_name(), "N");
    assert_eq!(labels[1].atoms(), ids(topology, &[1, 4]));
    assert!(error
        .to_string()
        .starts_with("2 moving label(s) name more than one selected atom; chain A residue 1 atom CA names atoms "));

    // The reference side is checked too, whether or not the moving side names
    // the label.
    let error =
        AtomCorrespondence::from_labels(&all(&reference), &all(&moving), LabelScheme::Author)
            .unwrap_err();
    assert!(matches!(
        error,
        AtomCorrespondenceError::AmbiguousLabels {
            side: CorrespondenceSide::Reference,
            ..
        }
    ));

    // Narrowing the selection resolves it.
    let first = AtomSelection::from_indices(
        moving.shared_topology(),
        [0, 1, 2].map(TopologyAtomIndex::new),
    )
    .unwrap();
    let paired =
        AtomCorrespondence::from_labels(&first, &all(&reference), LabelScheme::Author).unwrap();
    assert_eq!(paired_indices(&paired), [(0, 0)]);
}

fn read_mmcif(rows: &str) -> Model {
    let text = format!(
        "\
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
{rows}"
    );
    mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .into_model()
}

#[test]
fn label_scheme_reads_mmcif_label_identifiers() {
    // Author and label identifiers disagree: the moving structure numbers its
    // author residues from 101 and names its chain X.
    let moving = read_mmcif(
        "\
ATOM 1 N N GLY A 1 1 0.000 0.000 0.000 101 X N 1
ATOM 2 C CA GLY A 1 1 1.400 0.000 0.000 101 X CA 1
HETATM 3 O O HOH B 2 . 5.000 0.000 0.000 201 X O 1
",
    );
    let reference = read_mmcif(
        "\
ATOM 1 N N GLY A 1 1 0.000 0.000 0.000 1 A N 1
ATOM 2 C CA GLY A 1 1 1.400 0.000 0.000 1 A CA 1
HETATM 3 O O HOH B 2 . 5.000 0.000 0.000 201 X O 1
",
    );
    let by_label =
        AtomCorrespondence::from_labels(&all(&moving), &all(&reference), LabelScheme::Label)
            .unwrap();
    // Waters have no label sequence number, so they are unlabeled.
    assert_eq!(paired_indices(&by_label), [(0, 0), (1, 1)]);
    assert_eq!(unpaired_indices(&by_label, CorrespondenceSide::Moving), [2]);
    assert_eq!(
        unpaired_indices(&by_label, CorrespondenceSide::Reference),
        [2]
    );

    let by_author =
        AtomCorrespondence::from_labels(&all(&moving), &all(&reference), LabelScheme::Author)
            .unwrap();
    assert_eq!(paired_indices(&by_author), [(2, 2)]);
    assert_eq!(
        unpaired_indices(&by_author, CorrespondenceSide::Moving),
        [0, 1]
    );
}

#[test]
fn atoms_without_hierarchy_are_unpaired() {
    let moving = smiles::to_topology("CCO").unwrap();
    let reference = smiles::to_topology("CCO").unwrap();
    let paired = AtomCorrespondence::from_labels(
        &AtomSelection::all(moving),
        &AtomSelection::all(reference),
        LabelScheme::default(),
    )
    .unwrap();
    assert!(paired.correspondence().is_empty());
    assert_eq!(paired.unpaired_moving().len(), 3);
    assert_eq!(paired.unpaired_reference().len(), 3);
}
