//! Chemical component templates in the model of the wwPDB Chemical Component
//! Dictionary (CCD): residue and ligand definitions keyed by component
//! identifier, with named atoms, elements, formal charges, and bonds.
//!
//! Structural formats name each residue's component but rarely carry its bonds.
//! Templates are the interpretation context that supplies them.

mod bundled;
mod ccd;

use std::collections::BTreeMap;
use std::fmt;
use std::sync::{Arc, OnceLock};

use crate::core::Element;
use crate::diagnostics::{DiagnosticKind, SourceLocation};

use super::MmcifDocument;

/// Polymer linkage a component takes part in, from its dictionary type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ComponentLinking {
    /// Peptide-linking residues and termini (`L-PEPTIDE LINKING`, ...).
    Peptide,
    /// DNA and RNA linking residues and termini.
    Nucleotide,
    /// Saccharides.
    Saccharide,
    /// Non-polymer and other components.
    None,
}

impl ComponentLinking {
    fn from_type(kind: &str) -> Self {
        let kind = kind.to_ascii_uppercase();
        if kind.contains("PEPTIDE")
            && (kind.contains("LINKING") || kind.contains("TERMINUS"))
            && !kind.contains("PEPTIDE-LIKE")
        {
            Self::Peptide
        } else if kind.contains("DNA") || kind.contains("RNA") {
            Self::Nucleotide
        } else if kind.contains("SACCHARIDE") {
            Self::Saccharide
        } else {
            Self::None
        }
    }
}

/// One named atom of a component.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ComponentAtom {
    pub(crate) name: String,
    pub(crate) alternate_name: Option<String>,
    pub(crate) element: Element,
    pub(crate) isotope: Option<u16>,
    pub(crate) formal_charge: Option<i8>,
    pub(crate) aromatic: bool,
    pub(crate) leaving: bool,
}

impl ComponentAtom {
    pub fn name(&self) -> &str {
        &self.name
    }

    /// An older name for the same atom, such as `1HB` for `HB2`.
    pub fn alternate_name(&self) -> Option<&str> {
        self.alternate_name.as_deref()
    }

    pub const fn element(&self) -> Element {
        self.element
    }

    /// Mass number of an isotope-specific element symbol, such as `D`.
    pub const fn isotope(&self) -> Option<u16> {
        self.isotope
    }

    /// The template formal charge, or `None` where the dictionary states it as
    /// unknown.
    pub const fn formal_charge(&self) -> Option<i8> {
        self.formal_charge
    }

    pub const fn aromatic(&self) -> bool {
        self.aromatic
    }

    /// Whether the atom leaves when the component forms a polymer or other
    /// covalent link, such as the carboxyl `OXT` of an amino acid.
    pub const fn leaving(&self) -> bool {
        self.leaving
    }

    fn is_named(&self, name: &str) -> bool {
        self.name == name || self.alternate_name.as_deref() == Some(name)
    }
}

/// Bond order as the dictionary states it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ComponentBondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
    Delocalized,
    Pi,
    Polymeric,
}

/// One bond between two atoms of a component.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComponentBond {
    pub(crate) atoms: [usize; 2],
    pub(crate) order: ComponentBondOrder,
    pub(crate) aromatic: bool,
    pub(crate) metal_coordination: bool,
}

impl ComponentBond {
    /// Endpoints as indices into [`Component::atoms`].
    pub const fn atoms(&self) -> [usize; 2] {
        self.atoms
    }

    pub const fn order(&self) -> ComponentBondOrder {
        self.order
    }

    pub const fn aromatic(&self) -> bool {
        self.aromatic
    }

    /// Whether the dictionary flags the bond as metal coordination. Such a
    /// bond does not count toward the valence of its ligand atom.
    pub const fn metal_coordination(&self) -> bool {
        self.metal_coordination
    }
}

/// One component template: a dictionary entry, or a protonation variant
/// derived from one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Component {
    pub(crate) id: String,
    pub(crate) name: String,
    pub(crate) kind: String,
    pub(crate) variant: Option<String>,
    pub(crate) atoms: Vec<ComponentAtom>,
    pub(crate) bonds: Vec<ComponentBond>,
}

impl Component {
    /// The component identifier, such as `ALA` or `HEM`.
    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    /// The dictionary type, such as `L-PEPTIDE LINKING` or `NON-POLYMER`.
    pub fn kind(&self) -> &str {
        &self.kind
    }

    pub fn linking(&self) -> ComponentLinking {
        ComponentLinking::from_type(&self.kind)
    }

    /// The name of a protonation variant, such as `NE2-protonated` for a
    /// histidine tautomer; `None` for a dictionary entry.
    pub fn variant(&self) -> Option<&str> {
        self.variant.as_deref()
    }

    pub fn atoms(&self) -> &[ComponentAtom] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[ComponentBond] {
        &self.bonds
    }

    /// Finds an atom by its name or alternate name.
    pub fn atom(&self, name: &str) -> Option<(usize, &ComponentAtom)> {
        self.atoms
            .iter()
            .enumerate()
            .find(|(_, atom)| atom.name == name)
            .or_else(|| {
                self.atoms
                    .iter()
                    .enumerate()
                    .find(|(_, atom)| atom.is_named(name))
            })
    }

    /// Reads one dictionary entry from a CIF data block.
    pub fn from_mmcif_block(block: &super::MmcifBlock) -> Result<Self, ComponentError> {
        ccd::component_from_block(block)
    }
}

/// Failure to read a component definition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ComponentError {
    kind: DiagnosticKind,
    location: Option<SourceLocation>,
    message: String,
}

impl ComponentError {
    pub(crate) fn new(
        kind: DiagnosticKind,
        location: Option<SourceLocation>,
        message: impl Into<String>,
    ) -> Self {
        Self {
            kind,
            location,
            message: message.into(),
        }
    }

    pub const fn kind(&self) -> DiagnosticKind {
        self.kind
    }

    pub const fn location(&self) -> Option<SourceLocation> {
        self.location
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for ComponentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.location {
            Some(location) => write!(f, "component error at {location}: {}", self.message),
            None => write!(f, "component error: {}", self.message),
        }
    }
}

impl std::error::Error for ComponentError {}

#[derive(Debug, Clone)]
struct ComponentEntry {
    component: Arc<Component>,
    variants: Vec<Arc<Component>>,
}

/// Component templates keyed by identifier, with any protonation variants.
///
/// Libraries are cheap to clone; components are shared.
#[derive(Debug, Clone, Default)]
pub struct ComponentLibrary {
    entries: BTreeMap<String, ComponentEntry>,
}

impl ComponentLibrary {
    pub fn new() -> Self {
        Self::default()
    }

    /// The templates Kekule bundles: standard amino acids and nucleotides,
    /// water, common ions, cofactors, and crystallization additives, taken
    /// unmodified from the Chemical Component Dictionary, plus the neutral
    /// histidine tautomers derived from `HIS`.
    pub fn bundled() -> &'static Self {
        static BUNDLED: OnceLock<ComponentLibrary> = OnceLock::new();
        BUNDLED.get_or_init(bundled::library)
    }

    /// Reads every data block that defines component atoms, such as a CCD
    /// entry file or the complete `components.cif`.
    pub fn from_mmcif(document: &MmcifDocument) -> Result<Self, ComponentError> {
        let mut library = Self::new();
        for block in document.blocks() {
            if block.has_category("_chem_comp_atom") {
                library.insert(Component::from_mmcif_block(block)?);
            }
        }
        Ok(library)
    }

    /// Adds or replaces a component. Replacing one discards its variants,
    /// which were derived from the previous definition.
    pub fn insert(&mut self, component: Component) -> Option<Arc<Component>> {
        let id = component.id.to_ascii_uppercase();
        self.entries
            .insert(
                id,
                ComponentEntry {
                    component: Arc::new(component),
                    variants: Vec::new(),
                },
            )
            .map(|entry| entry.component)
    }

    /// Adds a protonation variant of a component already in the library. The
    /// variant must name the component and have a variant name.
    pub fn insert_variant(&mut self, variant: Component) -> Result<(), ComponentError> {
        let Some(name) = variant.variant.clone() else {
            return Err(ComponentError::new(
                DiagnosticKind::InvalidValue,
                None,
                format!("variant of `{}` has no variant name", variant.id),
            ));
        };
        let entry = self
            .entries
            .get_mut(&variant.id.to_ascii_uppercase())
            .ok_or_else(|| {
                ComponentError::new(
                    DiagnosticKind::InvalidReference,
                    None,
                    format!("variant `{name}` names unknown component `{}`", variant.id),
                )
            })?;
        entry
            .variants
            .retain(|existing| existing.variant.as_deref() != Some(name.as_str()));
        entry.variants.push(Arc::new(variant));
        Ok(())
    }

    /// Adds every component of `other`, with its variants, replacing
    /// components of the same identifier.
    pub fn extend(&mut self, other: &Self) {
        for (id, entry) in &other.entries {
            self.entries.insert(id.clone(), entry.clone());
        }
    }

    /// The component with identifier `id`, compared case-insensitively.
    pub fn get(&self, id: &str) -> Option<&Component> {
        self.entry(id).map(|entry| entry.component.as_ref())
    }

    /// Protonation variants of the component with identifier `id`.
    pub fn variants(&self, id: &str) -> impl ExactSizeIterator<Item = &Component> {
        self.entry(id)
            .map_or(&[][..], |entry| entry.variants.as_slice())
            .iter()
            .map(Arc::as_ref)
    }

    /// Components in identifier order.
    pub fn components(&self) -> impl ExactSizeIterator<Item = &Component> {
        self.entries.values().map(|entry| entry.component.as_ref())
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    fn entry(&self, id: &str) -> Option<&ComponentEntry> {
        self.entries.get(&id.to_ascii_uppercase())
    }
}
