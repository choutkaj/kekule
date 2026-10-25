//! Atom correspondence of two topologies by their hierarchy labels.

use std::collections::{BTreeMap, HashMap};
use std::fmt;

use super::{AtomCorrespondence, AtomCorrespondenceError, CorrespondenceSide};
use crate::topology::{AtomSelection, InstanceAtomId, Topology, TopologyAtomIndex};

/// Which hierarchy identifiers label an atom when pairing by labels.
///
/// The two schemes are the two namespaces mmCIF `atom_site` carries and kekule's
/// hierarchy keeps. Each reads exactly what the selection language's terms of
/// the same namespace read, so a label pairs the atoms those terms select.
///
/// This enum is closed: kekule's hierarchy has exactly these two namespaces.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub enum LabelScheme {
    /// Author identifiers, as PDB files and the selection language's `chain`,
    /// `resid`, and `name` terms read them. These are the labels that stay the
    /// same between depositions and models of one protein.
    ///
    /// - chain: the author chain ID, else the label chain ID;
    /// - residue number: the author sequence ID read as an integer, with the
    ///   residue's insertion code; a residue without an author sequence ID
    ///   takes its label sequence number and no insertion code;
    /// - atom name: the author atom ID, else the label atom ID.
    #[default]
    Author,
    /// mmCIF label identifiers, as the selection language's `labelchain`,
    /// `labelresid`, and `labelname` terms read them:
    /// `label_asym_id`, `label_seq_id`, and `label_atom_id`. Label sequence
    /// numbers have no insertion codes. Residues outside polymers, such as
    /// waters and ligands, usually have no label sequence number, so their
    /// atoms are unlabeled under this scheme.
    Label,
}

impl fmt::Display for LabelScheme {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Author => "author",
            Self::Label => "label",
        })
    }
}

/// The hierarchy label of one atom under a [`LabelScheme`]: chain, residue
/// number, insertion code, and atom name.
///
/// Residue names are not part of a label, so a mutated residue still pairs
/// the atoms both residues name, such as its backbone.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AtomLabel {
    chain: String,
    residue_number: i64,
    insertion_code: Option<String>,
    atom_name: String,
}

impl AtomLabel {
    /// The chain identifier.
    pub fn chain(&self) -> &str {
        &self.chain
    }

    /// The residue number.
    pub const fn residue_number(&self) -> i64 {
        self.residue_number
    }

    /// The insertion code, or `None` when the residue has none.
    pub fn insertion_code(&self) -> Option<&str> {
        self.insertion_code.as_deref()
    }

    /// The atom name.
    pub fn atom_name(&self) -> &str {
        &self.atom_name
    }
}

impl fmt::Display for AtomLabel {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "chain {} residue {}{} atom {}",
            self.chain,
            self.residue_number,
            self.insertion_code.as_deref().unwrap_or(""),
            self.atom_name,
        )
    }
}

/// One label that names more than one selected atom on one side.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AmbiguousLabel {
    label: AtomLabel,
    atoms: Vec<InstanceAtomId>,
}

impl AmbiguousLabel {
    /// The label shared by the atoms.
    pub fn label(&self) -> &AtomLabel {
        &self.label
    }

    /// The selected atoms the label names, two or more, in dense order.
    pub fn atoms(&self) -> &[InstanceAtomId] {
        &self.atoms
    }
}

impl fmt::Display for AmbiguousLabel {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{} names atoms ", self.label)?;
        for (position, atom) in self.atoms.iter().enumerate() {
            if position > 0 {
                formatter.write_str(", ")?;
            }
            write!(formatter, "{atom}")?;
        }
        Ok(())
    }
}

/// An [`AtomCorrespondence`] made by pairing hierarchy labels, with the
/// selected atoms of each side that it leaves unpaired.
///
/// See [`AtomCorrespondence::from_labels`].
#[derive(Debug, Clone)]
pub struct LabelCorrespondence {
    correspondence: AtomCorrespondence,
    unpaired_moving: AtomSelection,
    unpaired_reference: AtomSelection,
}

impl LabelCorrespondence {
    /// The atom pairs, in the moving selection's dense order.
    pub fn correspondence(&self) -> &AtomCorrespondence {
        &self.correspondence
    }

    /// The atom pairs, without the unpaired atoms.
    pub fn into_correspondence(self) -> AtomCorrespondence {
        self.correspondence
    }

    /// The selected atoms of one side left unpaired: those whose label the
    /// other side's selection does not name, and those without a label.
    pub fn unpaired(&self, side: CorrespondenceSide) -> &AtomSelection {
        match side {
            CorrespondenceSide::Moving => &self.unpaired_moving,
            CorrespondenceSide::Reference => &self.unpaired_reference,
        }
    }

    /// The selected moving atoms left unpaired.
    pub fn unpaired_moving(&self) -> &AtomSelection {
        &self.unpaired_moving
    }

    /// The selected reference atoms left unpaired.
    pub fn unpaired_reference(&self) -> &AtomSelection {
        &self.unpaired_reference
    }
}

impl AtomCorrespondence {
    /// Pairs the atoms of two selections whose hierarchy labels are equal.
    ///
    /// A label is an atom's chain, residue number, insertion code, and atom
    /// name, read under `scheme` ([`LabelScheme`] says which identifiers). Each
    /// selection's own topology snapshot is that side's topology, so the two
    /// sides may be independent structures: two crystal structures of one
    /// protein, or a model and its template. Pairs follow the moving
    /// selection's dense order. Elements, residue names, and bonding are not
    /// compared.
    ///
    /// Atoms without a complete label (no atom site, no atom name, or no
    /// residue number under `scheme`) and atoms whose label the other side's
    /// selection does not name stay unpaired; [`LabelCorrespondence::unpaired`]
    /// reports them. Labels are not trimmed or case-folded, except that a
    /// blank insertion code counts as none. Calculations enforce their own
    /// minimum number of pairs.
    ///
    /// # Errors
    ///
    /// [`AtomCorrespondenceError::AmbiguousLabels`] when a label names more
    /// than one selected atom on a side, for example alternate locations kept
    /// as separate atoms or repeated atom names. It names every such label of
    /// the first side, moving then reference, that has one, whether or not the
    /// other side names the label. Narrow the selection to resolve it.
    ///
    /// # Examples
    ///
    /// ```
    /// use kekule::alignment::{AtomCorrespondence, LabelScheme};
    /// use kekule::pdb;
    /// use kekule::topology::AtomSelection;
    ///
    /// let read = |text: &str| -> Result<_, Box<dyn std::error::Error>> {
    ///     Ok(pdb::parse_str(text, Default::default())?
    ///         .interpret(Default::default())?
    ///         .into_model())
    /// };
    /// // A model and a template of one water, the template with a second hydrogen.
    /// let model = read("\
    /// HETATM    1  O   HOH A   1       0.100   0.000   0.000  1.00  0.00           O
    /// HETATM    2  H1  HOH A   1       1.060   0.000   0.000  1.00  0.00           H
    /// END
    /// ")?;
    /// let template = read("\
    /// HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
    /// HETATM    2  H1  HOH A   1       0.960   0.000   0.000  1.00  0.00           H
    /// HETATM    3  H2  HOH A   1      -0.240   0.930   0.000  1.00  0.00           H
    /// END
    /// ")?;
    ///
    /// let paired = AtomCorrespondence::from_labels(
    ///     &AtomSelection::all(model.shared_topology()),
    ///     &AtomSelection::all(template.shared_topology()),
    ///     LabelScheme::Author,
    /// )?;
    /// assert_eq!(paired.correspondence().len(), 2);
    /// assert!(paired.unpaired_moving().is_empty());
    /// assert_eq!(paired.unpaired_reference().len(), 1);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn from_labels(
        moving: &AtomSelection,
        reference: &AtomSelection,
        scheme: LabelScheme,
    ) -> Result<LabelCorrespondence, AtomCorrespondenceError> {
        let moving_side = Side::new(moving, scheme, CorrespondenceSide::Moving)?;
        let reference_side = Side::new(reference, scheme, CorrespondenceSide::Reference)?;
        let mut reference_paired = vec![false; reference.len()];
        let mut pairs = Vec::new();
        let mut unpaired_moving = Vec::new();
        for (&index, key) in moving.indices().iter().zip(&moving_side.keys) {
            match key
                .as_ref()
                .and_then(|key| reference_side.positions.get(key))
            {
                Some(&position) => {
                    reference_paired[position] = true;
                    pairs.push((index, reference.indices()[position]));
                }
                None => unpaired_moving.push(index),
            }
        }
        let unpaired_reference = reference
            .indices()
            .iter()
            .zip(&reference_paired)
            .filter(|(_, &paired)| !paired)
            .map(|(&index, _)| index)
            .collect();
        let unpaired = |selection: &AtomSelection, indices| {
            AtomSelection::from_sorted_indices(selection.shared_topology(), indices)
                .expect("a subsequence of a selection's sorted indices")
        };
        Ok(LabelCorrespondence {
            correspondence: AtomCorrespondence {
                moving: moving.shared_topology(),
                reference: reference.shared_topology(),
                pairs,
            },
            unpaired_moving: unpaired(moving, unpaired_moving),
            unpaired_reference: unpaired(reference, unpaired_reference),
        })
    }
}

/// An atom label borrowed from its topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Key<'a> {
    chain: &'a str,
    number: i64,
    insertion: Option<&'a str>,
    name: &'a str,
}

impl<'a> Key<'a> {
    /// The label of the atom at `index`, or `None` when it has no complete one.
    fn of(topology: &'a Topology, index: TopologyAtomIndex, scheme: LabelScheme) -> Option<Self> {
        let site = topology.atom_site_at(index)?;
        let residue = site.residue();
        let chain = residue.chain();
        let metadata = site.metadata();
        let label_number = || residue.label_seq_id().map(i64::from);
        let (chain, number, insertion, name) = match scheme {
            LabelScheme::Author => {
                let (number, insertion) = match residue.author_seq_id() {
                    Some(text) => (
                        text.trim().parse().ok()?,
                        residue
                            .insertion_code()
                            .filter(|code| !code.trim().is_empty()),
                    ),
                    None => (label_number()?, None),
                };
                let name = metadata.auth_atom_id().or(metadata.label_atom_id())?;
                (
                    chain.author_id().unwrap_or(chain.label_id()),
                    number,
                    insertion,
                    name,
                )
            }
            LabelScheme::Label => (
                chain.label_id(),
                label_number()?,
                None,
                metadata.label_atom_id()?,
            ),
        };
        Some(Self {
            chain,
            number,
            insertion,
            name,
        })
    }

    fn to_label(self) -> AtomLabel {
        AtomLabel {
            chain: self.chain.to_owned(),
            residue_number: self.number,
            insertion_code: self.insertion.map(str::to_owned),
            atom_name: self.name.to_owned(),
        }
    }
}

/// The labels of one side's selected atoms.
struct Side<'a> {
    /// Each selected atom's label, in selection order.
    keys: Vec<Option<Key<'a>>>,
    /// The selection position of the one atom each label names.
    positions: HashMap<Key<'a>, usize>,
}

impl<'a> Side<'a> {
    fn new(
        selection: &'a AtomSelection,
        scheme: LabelScheme,
        side: CorrespondenceSide,
    ) -> Result<Self, AtomCorrespondenceError> {
        let topology = selection.topology();
        let keys: Vec<_> = selection
            .indices()
            .iter()
            .map(|&index| Key::of(topology, index, scheme))
            .collect();
        let mut positions = HashMap::with_capacity(keys.len());
        // Repeated labels by the position of their first atom, with the
        // positions of the atoms after it.
        let mut repeated: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (position, key) in keys.iter().enumerate() {
            if let Some(key) = key {
                let first = *positions.entry(*key).or_insert(position);
                if first != position {
                    repeated.entry(first).or_default().push(position);
                }
            }
        }
        if repeated.is_empty() {
            return Ok(Self { keys, positions });
        }
        let atom_id = |position: usize| {
            topology
                .atom_id(selection.indices()[position])
                .expect("selected atom of its own topology")
        };
        let labels = repeated
            .into_iter()
            .map(|(first, others)| AmbiguousLabel {
                label: keys[first].expect("labeled atom").to_label(),
                atoms: std::iter::once(first).chain(others).map(atom_id).collect(),
            })
            .collect();
        Err(AtomCorrespondenceError::AmbiguousLabels { side, labels })
    }
}

impl From<LabelCorrespondence> for AtomCorrespondence {
    fn from(paired: LabelCorrespondence) -> Self {
        paired.correspondence
    }
}
