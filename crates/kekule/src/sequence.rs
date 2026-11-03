//! One-letter codes of residues and the sequences of polymer chains.
//!
//! A residue's code comes from its chemical component, the component its
//! hierarchy names (`label_comp_id`, else `auth_comp_id`, else its name):
//!
//! 1. the code a component library states for it
//!    ([`SequenceOptions::components`], then the bundled templates), else the
//!    Chemical Component Dictionary's, which Kekule bundles for every
//!    component that links into a polypeptide or a polynucleotide;
//! 2. without one, the codes of the components it derives from
//!    (`mon_nstd_parent_comp_id`), so a modified residue takes its parent's
//!    code;
//! 3. without either, the explicit unknown: `X` for a peptide-linking
//!    component or an amino acid, `N` for a nucleotide.
//!
//! Only monomers have codes: residues whose components link into a
//! polypeptide or a polynucleotide, and amino-acid and nucleotide residues.
//! Other residues, such as water, ions and ligands, have none, even a ligand
//! such as `GTP` whose dictionary entry states one. Most codes
//! are one letter; a few components that stand for several monomers, such as
//! the chromophore of green fluorescent protein (`CRO`, `TYG`), have several.
//!
//! A chain's sequence lists its residues that have codes, in chain order, and
//! joins their codes. It is the sequence of the residues the topology holds;
//! the sequence a file states for a polymer entity, with the residues it does
//! not observe, comes with the mmCIF interpretation
//! ([`StructureInterpretationReport::polymer_entities`](crate::macromolecular::StructureInterpretationReport::polymer_entities)).
//!
//! ```
//! use kekule::sequence::{self, CodeOrigin, SequenceOptions};
//!
//! let options = SequenceOptions::default();
//! let threonine = sequence::component_code("TPO", &options).expect("phosphothreonine");
//! assert_eq!(threonine.as_str(), "T");
//! assert_eq!(threonine.origin(), CodeOrigin::Component);
//! assert!(sequence::component_code("HOH", &options).is_none());
//! ```

use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::sync::Arc;

use crate::io::{dictionary_code, ComponentLibrary, ComponentLinking};
use crate::topology::{ChainId, ResidueClass, ResidueId, ResidueView, Topology};

/// Where a residue's one-letter code comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum CodeOrigin {
    /// The code its component states.
    Component,
    /// The codes of the components it derives from, its parents.
    Parent,
    /// The canonical sequence of its polymer entity, as the file states it.
    /// Only monomers of an interpreted entity have it.
    Sequence,
    /// No code is known: `X` for a peptide, `N` for a nucleic acid.
    Unknown,
}

/// The one-letter code of a residue or monomer, usually one letter.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ResidueCode {
    code: Box<str>,
    origin: CodeOrigin,
}

impl ResidueCode {
    pub(crate) fn new(code: impl Into<Box<str>>, origin: CodeOrigin) -> Self {
        Self {
            code: code.into(),
            origin,
        }
    }

    /// The code, such as `A`, or `TYG` for a chromophore standing for three
    /// monomers.
    pub fn as_str(&self) -> &str {
        &self.code
    }

    /// Where the code comes from.
    pub const fn origin(&self) -> CodeOrigin {
        self.origin
    }

    /// Whether the code is the explicit unknown, `X` or `N`.
    pub fn is_unknown(&self) -> bool {
        self.origin == CodeOrigin::Unknown
    }
}

impl fmt::Display for ResidueCode {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.code)
    }
}

/// The component library that supplies codes.
#[derive(Debug, Clone, Default)]
pub struct SequenceOptions {
    components: Option<Arc<ComponentLibrary>>,
}

impl SequenceOptions {
    /// The caller's component library, consulted before the bundled templates
    /// and the dictionary's codes; `None` by default.
    pub fn components(&self) -> Option<&Arc<ComponentLibrary>> {
        self.components.as_ref()
    }

    /// Sets the caller's component library, such as
    /// [`ComponentLibrary::amber`] for structures with AMBER residue names.
    #[must_use]
    pub fn with_components(mut self, components: Option<Arc<ComponentLibrary>>) -> Self {
        self.components = components;
        self
    }
}

/// The code of the component `id` when it links into a polypeptide or a
/// polynucleotide: its own, else its parents', else the unknown. `None` for
/// any other component, or one no library knows.
pub fn component_code(id: &str, options: &SequenceOptions) -> Option<ResidueCode> {
    CodeSources::new(None, options.components.as_deref()).code(id, None)
}

/// The code of a residue, from its component; an amino acid or nucleotide
/// whose component has no code takes the unknown, `X` or `N`. `None` for
/// other residues.
pub fn residue_code(residue: ResidueView<'_>, options: &SequenceOptions) -> Option<ResidueCode> {
    CodeSources::new(None, options.components.as_deref()).code(
        residue_component(residue),
        unknown_of_class(residue.class()),
    )
}

/// The residues of one chain that have codes, in chain order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainSequence {
    chain: ChainId,
    monomers: Vec<ChainMonomer>,
    sequence: String,
}

impl ChainSequence {
    /// The chain.
    pub const fn chain(&self) -> ChainId {
        self.chain
    }

    /// The residues that have codes, in chain order.
    pub fn monomers(&self) -> &[ChainMonomer] {
        &self.monomers
    }

    /// The monomers' codes joined, such as `MKV…`.
    pub fn sequence(&self) -> &str {
        &self.sequence
    }
}

/// One residue of a [`ChainSequence`] and its code.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainMonomer {
    residue: ResidueId,
    code: ResidueCode,
}

impl ChainMonomer {
    /// The residue.
    pub const fn residue(&self) -> ResidueId {
        self.residue
    }

    /// Its code.
    pub const fn code(&self) -> &ResidueCode {
        &self.code
    }
}

/// The sequence of every chain with at least one residue that has a code, in
/// hierarchy order. Water, ions and ligands a chain holds are left out.
pub fn chain_sequences(topology: &Topology, options: &SequenceOptions) -> Vec<ChainSequence> {
    let sources = CodeSources::new(None, options.components.as_deref());
    let mut known = HashMap::<(&str, Option<&'static str>), Option<ResidueCode>>::new();
    let mut sequences = Vec::new();
    for chain in topology.chains() {
        let mut monomers = Vec::new();
        let mut sequence = String::new();
        for residue in chain.residues() {
            let unknown = unknown_of_class(residue.class());
            let component = residue_component(residue);
            let code = known
                .entry((component, unknown))
                .or_insert_with(|| sources.code(component, unknown));
            if let Some(code) = code {
                sequence.push_str(code.as_str());
                monomers.push(ChainMonomer {
                    residue: residue.id(),
                    code: code.clone(),
                });
            }
        }
        if !monomers.is_empty() {
            sequences.push(ChainSequence {
                chain: chain.id(),
                monomers,
                sequence,
            });
        }
    }
    sequences
}

fn residue_component<'a>(residue: ResidueView<'a>) -> &'a str {
    residue
        .label_comp_id()
        .or(residue.author_comp_id())
        .unwrap_or_else(|| residue.name())
}

const fn unknown_of_class(class: ResidueClass) -> Option<&'static str> {
    match class {
        ResidueClass::AminoAcid => Some("X"),
        ResidueClass::DnaNucleotide | ResidueClass::RnaNucleotide => Some("N"),
        _ => None,
    }
}

/// The unknown code of a component that links into a polymer.
pub(crate) const fn unknown_of_linking(linking: ComponentLinking) -> Option<&'static str> {
    match linking {
        ComponentLinking::Peptide => Some("X"),
        ComponentLinking::Nucleotide => Some("N"),
        _ => None,
    }
}

/// A component as a source file's own `_chem_comp` category states it.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct StatedComponent {
    /// The linkage of its `type`, when stated.
    pub(crate) linking: Option<ComponentLinking>,
    pub(crate) code: Option<String>,
    pub(crate) parents: Vec<String>,
}

/// What one source states about a component.
#[derive(Default)]
struct Facts<'a> {
    code: Option<&'a str>,
    parents: Vec<&'a str>,
    linking: Option<ComponentLinking>,
}

/// The sources of codes in precedence order: a file's own components, the
/// caller's library, the bundled templates, and the dictionary's codes.
pub(crate) struct CodeSources<'a> {
    stated: Option<&'a BTreeMap<String, StatedComponent>>,
    library: Option<&'a ComponentLibrary>,
}

/// How many generations of parents a code is looked up through.
const PARENT_DEPTH: usize = 4;

impl<'a> CodeSources<'a> {
    /// `stated` is keyed by upper-case component identifier.
    pub(crate) const fn new(
        stated: Option<&'a BTreeMap<String, StatedComponent>>,
        library: Option<&'a ComponentLibrary>,
    ) -> Self {
        Self { stated, library }
    }

    /// The code of component `id`, or `unknown` when it has none and its
    /// linkage does not say which unknown applies. Only a monomer has a code:
    /// a component that links into a polypeptide or a polynucleotide, or one
    /// whose context supplies an unknown, such as an amino-acid residue or a
    /// monomer of a polymer entity. Ligands such as `GTP`, whose dictionary
    /// entries state codes, have none.
    pub(crate) fn code(&self, id: &str, unknown: Option<&'static str>) -> Option<ResidueCode> {
        let facts = self.facts(id);
        let linked = facts.linking.and_then(unknown_of_linking);
        if linked.is_none() && unknown.is_none() {
            return None;
        }
        if let Some(code) = facts.code {
            return Some(ResidueCode::new(code, CodeOrigin::Component));
        }
        if let Some(code) = self.parent_code(&facts.parents, PARENT_DEPTH) {
            return Some(ResidueCode::new(code, CodeOrigin::Parent));
        }
        linked
            .or(unknown)
            .map(|code| ResidueCode::new(code, CodeOrigin::Unknown))
    }

    /// The linkage the first source that knows `id` states.
    pub(crate) fn linking(&self, id: &str) -> Option<ComponentLinking> {
        self.facts(id).linking
    }

    fn parent_code(&self, parents: &[&str], depth: usize) -> Option<String> {
        if parents.is_empty() || depth == 0 {
            return None;
        }
        let mut code = String::new();
        for parent in parents {
            let facts = self.facts(parent);
            match facts.code {
                Some(own) => code.push_str(own),
                None => code.push_str(&self.parent_code(&facts.parents, depth - 1)?),
            }
        }
        Some(code)
    }

    /// Each fact from the first source that states it.
    fn facts(&self, id: &str) -> Facts<'a> {
        let mut facts = Facts::default();
        let mut add = |code: Option<&'a str>, parents: Vec<&'a str>, linking| {
            facts.code = facts.code.or(code);
            if facts.parents.is_empty() {
                facts.parents = parents;
            }
            facts.linking = facts.linking.or(linking);
        };
        let id = id.trim();
        if let Some(stated) = self
            .stated
            .and_then(|stated| stated.get(&id.to_ascii_uppercase()))
        {
            add(
                stated.code.as_deref(),
                stated.parents.iter().map(String::as_str).collect(),
                stated.linking,
            );
        }
        for library in self
            .library
            .into_iter()
            .chain([ComponentLibrary::bundled()])
        {
            if let Some(component) = library.get(id) {
                add(
                    component.one_letter_code(),
                    component.parent_ids().iter().map(String::as_str).collect(),
                    Some(component.linking()),
                );
            }
        }
        if let Some(entry) = dictionary_code(id) {
            add(entry.code, entry.parents().collect(), Some(entry.linking));
        }
        facts
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn code(id: &str) -> Option<(String, CodeOrigin)> {
        component_code(id, &SequenceOptions::default())
            .map(|code| (code.as_str().to_owned(), code.origin()))
    }

    #[test]
    fn standard_and_modified_components() {
        assert_eq!(code("ALA"), Some(("A".into(), CodeOrigin::Component)));
        assert_eq!(code("dc"), Some(("C".into(), CodeOrigin::Component)));
        assert_eq!(code("MSE"), Some(("M".into(), CodeOrigin::Component)));
        assert_eq!(code("CRO"), Some(("TYG".into(), CodeOrigin::Component)));
        // The dictionary states no code for 02O, a modified alanine.
        assert_eq!(code("02O"), Some(("A".into(), CodeOrigin::Parent)));
        // A dinucleotide with no code of its own takes both parents' codes.
        assert_eq!(code("05A"), Some(("TT".into(), CodeOrigin::Parent)));
        // 004 links into peptides but has neither a code nor a parent.
        assert_eq!(code("004"), Some(("X".into(), CodeOrigin::Unknown)));
        assert_eq!(code("HOH"), None);
        assert_eq!(code("NOT-A-COMPONENT"), None);
        // Nucleotide ligands state codes or parents but link into no polymer.
        assert_eq!(code("GTP"), None);
        assert_eq!(code("ADP"), None);
    }

    #[test]
    fn stated_components_come_first() {
        let stated = BTreeMap::from([
            (
                "ALA".to_owned(),
                StatedComponent {
                    linking: Some(ComponentLinking::Peptide),
                    ..StatedComponent::default()
                },
            ),
            (
                "XYZ".to_owned(),
                StatedComponent {
                    linking: Some(ComponentLinking::Nucleotide),
                    parents: vec!["DA".to_owned()],
                    ..StatedComponent::default()
                },
            ),
            (
                "QQQ".to_owned(),
                StatedComponent {
                    linking: Some(ComponentLinking::Nucleotide),
                    ..StatedComponent::default()
                },
            ),
            (
                "LOOP".to_owned(),
                StatedComponent {
                    linking: Some(ComponentLinking::Peptide),
                    parents: vec!["LOOP".to_owned()],
                    ..StatedComponent::default()
                },
            ),
        ]);
        let sources = CodeSources::new(Some(&stated), None);
        // A stated component without a code still takes the libraries' code.
        assert_eq!(
            sources.code("ALA", None).map(|code| code.code),
            Some("A".into())
        );
        let modified = sources.code("xyz", None).expect("stated parent");
        assert_eq!(
            (modified.as_str(), modified.origin()),
            ("A", CodeOrigin::Parent)
        );
        let unknown = sources.code("QQQ", Some("X")).expect("stated nucleotide");
        assert_eq!((unknown.as_str(), unknown.is_unknown()), ("N", true));
        assert_eq!(
            sources.code("LOOP", None).map(|code| code.code),
            Some("X".into())
        );
    }
}
