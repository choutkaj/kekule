//! Resolves the bonds a macromolecular source does not state, from ranked
//! evidence: the source's own component bond table or connection records,
//! component templates, polymer and branch linkage, and, for residues nothing
//! else covers, covalent radii.
use crate::diagnostics::{Cancellation, DiagnosticKind, FormatError};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use crate::algorithms::{
    implied_charges, incomplete_metal_valence, is_metal_donor, onium_hydrogens, ImpliedCharges,
};
use crate::components::{
    Component, ComponentBond as TemplateBond, ComponentBondOrder as TemplateOrder,
    ComponentLibrary, ComponentLinking,
};
use crate::io::localize_source_aromatic_bonds;
use crate::math::Point3;
use crate::molecule::{AtomId, BondOrder, Element, ImplicitHydrogens, MoleculeEditor};

use crate::io::macromolecular as raw;
use crate::io::mmcif::category::MmcifCategory;
use crate::io::MmcifBlock;

/// One staged atom with the source fields connectivity resolution reads.
pub(crate) struct StagedAtomProvenance<'a> {
    pub(crate) atom: AtomId,
    pub(crate) atom_name: &'a str,
    pub(crate) component_id: &'a str,
    pub(crate) asym_id: &'a str,
    pub(crate) entity_id: Option<&'a str>,
    pub(crate) label_sequence_id: Option<i32>,
    /// The residue's position along its polymer chain; residues adjacent in
    /// position are adjacent in sequence.
    pub(crate) polymer_position: Option<i32>,
    pub(crate) author_sequence_id: Option<&'a str>,
    pub(crate) insertion_code: Option<&'a str>,
    pub(crate) occurrence: Option<usize>,
    pub(crate) element: Element,
    pub(crate) point: Point3,
    pub(crate) line: usize,
    pub(crate) stated_charge: Option<i8>,
}

/// The evidence of each staged bond, keyed by its ordered atom pair. The first
/// evidence recorded for a pair is kept.
pub(crate) type BondEvidenceMap = BTreeMap<(AtomId, AtomId), raw::BondEvidence>;

/// A bond a connection record states between two atoms of one residue. It
/// is checked against the residue's template, and replaces distance
/// inference for residues without one.
#[derive(Debug, Clone)]
pub(crate) struct ListedBond {
    pub(crate) atoms: [AtomId; 2],
    pub(crate) order: BondOrder,
    pub(crate) order_source: OrderSource,
    pub(crate) evidence: raw::BondEvidence,
}

/// How a connection record gives the order of its bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum OrderSource {
    /// The record states the order, as PDB `SSBOND` records and mmCIF
    /// `_struct_conn` rows with `pdbx_value_order` do. A stated order
    /// outranks a template's, which the report notes.
    Stated,
    /// A listing convention gives the order: PDB `CONECT` records list a
    /// double bond twice, and a bond between a metal and a non-metal is
    /// coordination. A template's order outranks it.
    Convention,
    /// The record gives no order, and the bond is single.
    Unstated,
}

pub(crate) fn record_evidence(
    evidence: &mut BondEvidenceMap,
    a: AtomId,
    b: AtomId,
    source: raw::BondEvidence,
) {
    evidence.entry((a.min(b), a.max(b))).or_insert(source);
}

/// Resolves the bonds of one staged instance in place. Bonds already
/// present, from connection records, must agree with every other source.
#[allow(clippy::too_many_arguments)]
pub(crate) fn complete_editor_connectivity<'a>(
    catalog: &ConnectivityCatalog,
    editor: &mut MoleculeEditor,
    atoms: impl IntoIterator<Item = StagedAtomProvenance<'a>>,
    listed: &[ListedBond],
    model_id: &str,
    evidence: &mut BondEvidenceMap,
    issues: &mut Vec<raw::StructureInterpretIssue>,
    cancellation: Option<&Cancellation>,
) -> Result<(), FormatError> {
    let mut source_aromatic_bonds = BTreeSet::new();
    let mut resolver = Resolver {
        catalog,
        graph: &mut *editor,
        source_aromatic_bonds: &mut source_aromatic_bonds,
        evidence,
        issues,
        model_id,
        cancellation,
    };
    let residues = residue_atoms(atoms, cancellation)?;
    resolver.resolve(residues, listed)?;
    raw::poll(cancellation)?;
    localize_source_aromatic_bonds(editor, &source_aromatic_bonds).map_err(interpret_error)
}

/// The linkage between two residues adjacent in a polymer chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PolymerLink {
    /// The peptide bond from `C` to the next residue's `N`.
    Peptide,
    /// The phosphodiester bond from `O3'` (or `O3*`) to the next residue's
    /// `P`.
    Nucleotide,
    /// A peptide bond, or else a phosphodiester bond, by atom names.
    Either,
}

impl PolymerLink {
    /// The linkage of two adjacent residues without a declared polymer type:
    /// their templates name it. A residue without a polymer template defers
    /// to its neighbour, and two such residues link by atom names; bonding
    /// distance decides either way.
    pub(crate) fn undeclared(left: ComponentLinking, right: ComponentLinking) -> Option<Self> {
        let polymer = |linking| {
            matches!(
                linking,
                ComponentLinking::Peptide | ComponentLinking::Nucleotide
            )
            .then_some(linking)
        };
        match (polymer(left), polymer(right)) {
            (Some(first), Some(second)) if first != second => None,
            (Some(ComponentLinking::Peptide), _) | (_, Some(ComponentLinking::Peptide)) => {
                Some(Self::Peptide)
            }
            (Some(_), _) | (_, Some(_)) => Some(Self::Nucleotide),
            (None, None) => Some(Self::Either),
        }
    }

    /// The linking atoms of the earlier and the later residue, which `left`
    /// and `right` look up by atom name.
    pub(crate) fn atoms<T>(
        self,
        left: impl Fn(&str) -> Option<T>,
        right: impl Fn(&str) -> Option<T>,
    ) -> Option<(T, T)> {
        let peptide = || left("C").zip(right("N"));
        let nucleotide = || left("O3'").or_else(|| left("O3*")).zip(right("P"));
        match self {
            Self::Peptide => peptide(),
            Self::Nucleotide => nucleotide(),
            Self::Either => peptide().or_else(nucleotide),
        }
    }
}

/// What a reader makes of the standard link between two residues adjacent
/// along a chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum StandardLink<A> {
    /// The earlier and the later residue's linking atoms are bonded.
    Linked([A; 2]),
    /// A residue lacks the linking atom.
    Missing,
    /// A bond joins a linking atom to an atom of another residue than its
    /// partner, as in an isopeptide or an N-methylated backbone.
    Skipped,
    /// The linking atoms are beyond bonding distance: a chain break.
    TooFar,
}

/// The reader's rule for the standard `link` between two residues adjacent
/// along a chain, which readers and writers share, so that a writer leaves
/// out exactly the links a reader rebuilds. `left` and `right` find the
/// earlier and the later residue's atoms by name. `bonded_elsewhere` tells
/// whether a bond joins an atom to one of another residue than its partner,
/// coordination aside, and `within_distance` whether two atoms lie within
/// bonding distance.
pub(crate) fn standard_link<A: Copy, E>(
    link: PolymerLink,
    left: impl Fn(&str) -> Option<A>,
    right: impl Fn(&str) -> Option<A>,
    mut bonded_elsewhere: impl FnMut(A, A) -> Result<bool, E>,
    within_distance: impl FnOnce(A, A) -> bool,
) -> Result<StandardLink<A>, E> {
    let Some((left, right)) = link.atoms(left, right) else {
        return Ok(StandardLink::Missing);
    };
    if bonded_elsewhere(left, right)? || bonded_elsewhere(right, left)? {
        return Ok(StandardLink::Skipped);
    }
    Ok(if within_distance(left, right) {
        StandardLink::Linked([left, right])
    } else {
        StandardLink::TooFar
    })
}

/// Whether a bond between residues is coordination, which readers record
/// as connection metadata rather than bond: dative, between a metal and a
/// non-metal.
pub(crate) fn is_coordination(order: BondOrder, elements: [Element; 2]) -> bool {
    order == BondOrder::Dative && raw::is_metal(elements[0]) != raw::is_metal(elements[1])
}

/// The bond of `candidate`, a template or one of its variants, between the
/// atoms of a stated bond of a residue read with `template`.
fn candidate_bond<'c>(
    template: &Component,
    candidate: &'c Component,
    bond: &StatedBond,
) -> Option<&'c TemplateBond> {
    let [a, b] = bond
        .template_atoms
        .map(|index| candidate.atom(template.atoms()[index].name()));
    let (Some((a, _)), Some((b, _))) = (a, b) else {
        return None;
    };
    candidate.bonds().iter().find(|candidate| {
        let [x, y] = candidate.atoms();
        (x, y) == (a, b) || (x, y) == (b, a)
    })
}

/// Whether `candidate`, a template or one of its variants, gives a stated
/// bond its order.
fn candidate_gives(template: &Component, candidate: &Component, bond: &StatedBond) -> bool {
    candidate_bond(template, candidate, bond)
        .is_some_and(|template_bond| gives_order(template_bond, bond.order))
}

/// Whether a template bond gives a stated order: aromatic template bonds
/// accept either localized order, and coordination is dative.
fn gives_order(bond: &TemplateBond, order: BondOrder) -> bool {
    match SourceBondOrder::of_template(bond) {
        SourceBondOrder::Localized(template) => template == order,
        SourceBondOrder::Aromatic => matches!(order, BondOrder::Single | BondOrder::Double),
    }
}

/// The atoms polymer links join: peptide `N` and `C`, nucleotide `P` and
/// `O3'`, also written `O3*`.
const LINKING_ATOMS: [&str; 5] = ["N", "C", "P", "O3'", "O3*"];

/// Whether a residue has an atom a polymer link joins.
fn has_linking_atom(residue: &ResidueAtoms) -> bool {
    LINKING_ATOMS
        .iter()
        .any(|name| residue.atoms.contains_key(*name))
}

/// The bond a template gives its two atoms, or `None` for an aromatic or
/// delocalized bond, whose order localization chooses.
pub(crate) const fn template_bond_order(bond: &TemplateBond) -> Option<BondOrder> {
    match SourceBondOrder::of_template(bond) {
        SourceBondOrder::Localized(order) => Some(order),
        SourceBondOrder::Aromatic => None,
    }
}

/// A bond order as a source states it, before aromatic localization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SourceBondOrder {
    Localized(BondOrder),
    Aromatic,
}

impl SourceBondOrder {
    const fn staged_order(self) -> BondOrder {
        match self {
            Self::Localized(order) => order,
            Self::Aromatic => BondOrder::Single,
        }
    }

    /// Metal coordination within a component does not count toward the
    /// valence of the ligand atom, so it is a dative bond whatever its stated
    /// order.
    const fn of_template(bond: &TemplateBond) -> Self {
        if bond.metal_coordination() {
            return Self::Localized(BondOrder::Dative);
        }
        match bond.order() {
            TemplateOrder::Double | TemplateOrder::Pi => Self::Localized(BondOrder::Double),
            TemplateOrder::Triple => Self::Localized(BondOrder::Triple),
            TemplateOrder::Quadruple => Self::Localized(BondOrder::Quadruple),
            TemplateOrder::Aromatic | TemplateOrder::Delocalized => Self::Aromatic,
            _ => Self::Localized(BondOrder::Single),
        }
    }
}

#[derive(Debug, Clone)]
struct ComponentBond {
    atom_1: String,
    atom_2: String,
    order: SourceBondOrder,
    /// The row's `pdbx_metal_coordination_flag`, or `None` when the source
    /// does not state it, as wwPDB entry files do not.
    coordination: Option<bool>,
}

#[derive(Debug, Clone)]
struct BranchSchemeSite {
    entity_id: String,
    asym_id: String,
    number: i32,
    component_id: String,
    author_sequence_id: Option<String>,
}

#[derive(Debug, Clone)]
struct BranchLink {
    entity_id: String,
    number_1: i32,
    atom_1: String,
    number_2: i32,
    atom_2: String,
    order: SourceBondOrder,
}

/// Block-level connectivity evidence and the templates to resolve residues with.
#[derive(Debug, Default)]
pub(crate) struct ConnectivityCatalog {
    component_bonds: BTreeMap<String, Vec<ComponentBond>>,
    polymer_types: BTreeMap<String, String>,
    branch_sites: Vec<BranchSchemeSite>,
    branch_links: Vec<BranchLink>,
    components: Option<Arc<ComponentLibrary>>,
    infer_bonds_from_distance: bool,
}

impl ConnectivityCatalog {
    /// A catalog of templates alone, for sources without block-level evidence.
    pub(crate) fn new(options: &raw::ConnectivityOptions) -> Self {
        Self {
            components: options.components().cloned(),
            infer_bonds_from_distance: options.infer_bonds_from_distance(),
            ..Self::default()
        }
    }

    pub(crate) fn from_block(
        block: &MmcifBlock,
        options: &raw::ConnectivityOptions,
    ) -> Result<Self, FormatError> {
        let mut catalog = Self::new(options);
        if let Some(table) = block.category("_chem_comp_bond")? {
            let table = &table;
            for row in 0..table.row_count() {
                let comp_id = required(table, row, "_chem_comp_bond.comp_id")?;
                let atom_1 = required(table, row, "_chem_comp_bond.atom_id_1")?;
                let atom_2 = required(table, row, "_chem_comp_bond.atom_id_2")?;
                let coordination =
                    optional(table, row, "_chem_comp_bond.pdbx_metal_coordination_flag")
                        .map(is_yes);
                let order = component_bond_order(
                    optional(table, row, "_chem_comp_bond.value_order").unwrap_or("sing"),
                    table,
                    row,
                )?;
                catalog
                    .component_bonds
                    .entry(normalized(comp_id))
                    .or_default()
                    .push(ComponentBond {
                        atom_1: atom_1.to_owned(),
                        atom_2: atom_2.to_owned(),
                        order,
                        coordination,
                    });
            }
        }
        if let Some(table) = block.category("_entity_poly")? {
            let table = &table;
            for row in 0..table.row_count() {
                let entity = required(table, row, "_entity_poly.entity_id")?;
                let kind = required(table, row, "_entity_poly.type")?;
                catalog
                    .polymer_types
                    .insert(entity.to_owned(), kind.to_owned());
            }
        }
        if let Some(table) = block.category("_pdbx_branch_scheme")? {
            let table = &table;
            for row in 0..table.row_count() {
                let entity_id = required(table, row, "_pdbx_branch_scheme.entity_id")?;
                let asym_id = required(table, row, "_pdbx_branch_scheme.asym_id")?;
                let number = required_i32(table, row, "_pdbx_branch_scheme.num")?;
                let component_id = required(table, row, "_pdbx_branch_scheme.mon_id")?;
                let author_sequence_id = optional(table, row, "_pdbx_branch_scheme.pdb_seq_num")
                    .or_else(|| optional(table, row, "_pdbx_branch_scheme.auth_seq_num"))
                    .map(str::to_owned);
                catalog.branch_sites.push(BranchSchemeSite {
                    entity_id: entity_id.to_owned(),
                    asym_id: asym_id.to_owned(),
                    number,
                    component_id: component_id.to_owned(),
                    author_sequence_id,
                });
            }
        }
        if let Some(table) = block.category("_pdbx_entity_branch_link")? {
            let table = &table;
            for row in 0..table.row_count() {
                let entity_id = required(table, row, "_pdbx_entity_branch_link.entity_id")?;
                let number_1 = required_i32(
                    table,
                    row,
                    "_pdbx_entity_branch_link.entity_branch_list_num_1",
                )?;
                let atom_1 = required(table, row, "_pdbx_entity_branch_link.atom_id_1")?;
                let number_2 = required_i32(
                    table,
                    row,
                    "_pdbx_entity_branch_link.entity_branch_list_num_2",
                )?;
                let atom_2 = required(table, row, "_pdbx_entity_branch_link.atom_id_2")?;
                let order = component_bond_order(
                    optional(table, row, "_pdbx_entity_branch_link.value_order").unwrap_or("sing"),
                    table,
                    row,
                )?;
                catalog.branch_links.push(BranchLink {
                    entity_id: entity_id.to_owned(),
                    number_1,
                    atom_1: atom_1.to_owned(),
                    number_2,
                    atom_2: atom_2.to_owned(),
                    order,
                });
            }
        }
        Ok(catalog)
    }

    /// The template of a component and its protonation variants: the caller's
    /// library takes precedence over the bundled one.
    fn template(&self, component_id: &str) -> Option<(&Component, Vec<&Component>)> {
        let library = match self.components.as_deref() {
            Some(library) if library.get(component_id).is_some() => library,
            _ => ComponentLibrary::bundled(),
        };
        let component = library.get(component_id)?;
        Some((component, library.variants(component_id).collect()))
    }

    /// The template a residue is read with: the caller's of its name, else
    /// the default one for its atoms (see [`ComponentLibrary::default_template`]),
    /// with the force field that names it, if any.
    fn residue_template(&self, residue: &ResidueAtoms) -> Option<ResidueTemplate<'_>> {
        let component_id = normalized(&residue.key.component_id);
        if let Some(library) = self
            .components
            .as_deref()
            .filter(|library| library.get(&component_id).is_some())
        {
            let component = library.get(&component_id)?;
            return Some(ResidueTemplate {
                component,
                variants: library.variants(&component_id).collect(),
                force_field: None,
            });
        }
        let atoms = residue
            .sites
            .iter()
            .map(|site| (site.name.as_str(), site.element));
        let template = ComponentLibrary::default_template(&component_id, atoms)?;
        Some(ResidueTemplate {
            component: template.component,
            variants: template.library.variants(&component_id).collect(),
            force_field: template.force_field,
        })
    }

    /// Whether bonds no evidence states are inferred from distances.
    pub(crate) const fn infers_bonds_from_distance(&self) -> bool {
        self.infer_bonds_from_distance
    }

    fn linking(&self, component_id: &str) -> ComponentLinking {
        self.template(component_id)
            .map_or(ComponentLinking::None, |(component, _)| component.linking())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ResidueKey {
    asym_id: String,
    entity_id: Option<String>,
    label_sequence_id: Option<i32>,
    author_sequence_id: Option<String>,
    insertion_code: Option<String>,
    occurrence: Option<usize>,
    component_id: String,
    polymer_position: Option<i32>,
}

impl ResidueKey {
    /// Whether `atom` belongs to the residue this key names.
    fn matches(&self, atom: &StagedAtomProvenance<'_>) -> bool {
        self.label_sequence_id == atom.label_sequence_id
            && self.polymer_position == atom.polymer_position
            && self.occurrence == atom.occurrence
            && self.asym_id == atom.asym_id
            && self.component_id == atom.component_id
            && self.entity_id.as_deref() == atom.entity_id
            && self.author_sequence_id.as_deref() == atom.author_sequence_id
            && self.insertion_code.as_deref() == atom.insertion_code
    }

    fn id(&self, model_id: &str) -> raw::SourceResidueId {
        raw::SourceResidueId::new(
            model_id,
            &self.asym_id,
            self.label_sequence_id,
            self.author_sequence_id.as_deref(),
            self.insertion_code.as_deref(),
            self.occurrence,
            &self.component_id,
        )
    }
}

#[derive(Debug, Clone)]
struct ResidueSite {
    name: String,
    atom: AtomId,
    element: Element,
    point: Point3,
    line: usize,
    stated_charge: Option<i8>,
}

#[derive(Debug, Clone)]
struct ResidueAtoms {
    key: ResidueKey,
    atoms: BTreeMap<String, AtomId>,
    sites: Vec<ResidueSite>,
}

impl ResidueAtoms {
    fn site(&self, atom: AtomId) -> Option<&ResidueSite> {
        self.sites.iter().find(|site| site.atom == atom)
    }
}

/// The template a residue is read with, its protonation variants, and the
/// force field that names it, if it is not a dictionary template.
struct ResidueTemplate<'c> {
    component: &'c Component,
    variants: Vec<&'c Component>,
    force_field: Option<&'static str>,
}

/// A residue whose bonds came from a template; `matched` maps template atom
/// indices to staged atoms.
struct TemplatedResidue<'c> {
    key: ResidueKey,
    component: &'c Component,
    variants: Vec<&'c Component>,
    matched: BTreeMap<usize, AtomId>,
    /// Template bonds whose order a record states, checked against the
    /// template or variant the residue takes.
    stated: Vec<StatedBond>,
}

/// A template bond whose order a connection record states.
struct StatedBond {
    atoms: [AtomId; 2],
    /// The template atom indices.
    template_atoms: [usize; 2],
    names: [String; 2],
    order: BondOrder,
}

fn residue_atoms<'a>(
    atoms: impl IntoIterator<Item = StagedAtomProvenance<'a>>,
    cancellation: Option<&Cancellation>,
) -> Result<BTreeMap<ResidueKey, ResidueAtoms>, FormatError> {
    let mut residues = BTreeMap::new();
    // A residue's atoms usually follow each other, so the residue receiving
    // atoms is looked up once per run of them.
    let mut current: Option<ResidueAtoms> = None;
    for (step, atom) in atoms.into_iter().enumerate() {
        raw::poll_every(cancellation, step)?;
        let residue = match &mut current {
            Some(residue) if residue.key.matches(&atom) => residue,
            _ => {
                if let Some(residue) = current.take() {
                    residues.insert(residue.key.clone(), residue);
                }
                let key = ResidueKey {
                    asym_id: atom.asym_id.to_owned(),
                    entity_id: atom.entity_id.map(str::to_owned),
                    label_sequence_id: atom.label_sequence_id,
                    author_sequence_id: atom.author_sequence_id.map(str::to_owned),
                    insertion_code: atom.insertion_code.map(str::to_owned),
                    occurrence: atom.occurrence,
                    component_id: atom.component_id.to_owned(),
                    polymer_position: atom.polymer_position,
                };
                current.insert(residues.remove(&key).unwrap_or(ResidueAtoms {
                    key,
                    atoms: BTreeMap::new(),
                    sites: Vec::new(),
                }))
            }
        };
        residue.atoms.insert(atom.atom_name.to_owned(), atom.atom);
        residue.sites.push(ResidueSite {
            name: atom.atom_name.to_owned(),
            atom: atom.atom,
            element: atom.element,
            point: atom.point,
            line: atom.line,
            stated_charge: atom.stated_charge,
        });
    }
    if let Some(residue) = current {
        residues.insert(residue.key.clone(), residue);
    }
    raw::poll(cancellation)?;
    Ok(residues)
}

/// The distance, in ångström, that distance inference adds to the sum of two
/// covalent radii.
pub(crate) const COVALENT_BOND_TOLERANCE_ANGSTROM: f64 = 0.45;
/// The distance, in ångström, at or below which two atoms are too close to
/// bond.
pub(crate) const MIN_COVALENT_BOND_DISTANCE_ANGSTROM: f64 = 0.4;
/// The covalent radius, in ångström, assumed for an element without a
/// tabulated one.
pub(crate) const FALLBACK_COVALENT_RADIUS_ANGSTROM: f64 = 0.77;

/// The covalent radius distance inference uses for `element`, in ångström,
/// and whether it is the fallback for an element without a tabulated one.
pub(crate) fn covalent_bonding_radius(element: Element) -> (f64, bool) {
    match element.covalent_radius_angstrom() {
        Some(radius) => (radius, false),
        None => (FALLBACK_COVALENT_RADIUS_ANGSTROM, true),
    }
}

/// Whether two atoms lie within covalent bonding distance: the sum of their
/// covalent radii plus a tolerance, and not implausibly close.
pub(crate) fn within_bonding_distance(
    (left, left_point): (Element, Point3),
    (right, right_point): (Element, Point3),
) -> bool {
    let cutoff = covalent_bonding_radius(left).0
        + covalent_bonding_radius(right).0
        + COVALENT_BOND_TOLERANCE_ANGSTROM;
    let distance_squared = (left_point.x - right_point.x).powi(2)
        + (left_point.y - right_point.y).powi(2)
        + (left_point.z - right_point.z).powi(2);
    distance_squared > MIN_COVALENT_BOND_DISTANCE_ANGSTROM.powi(2)
        && distance_squared <= cutoff * cutoff
}

/// The coordination flag `component` gives the bond between the two named
/// atoms of a component-table row, if it lists that bond and its atoms have
/// the source atoms' elements; templates match by name and element.
fn template_coordination(
    component: &Component,
    bond: &ComponentBond,
    elements: [Element; 2],
) -> Option<bool> {
    let (left, left_atom) = component.atom(&bond.atom_1)?;
    let (right, right_atom) = component.atom(&bond.atom_2)?;
    if [left_atom.element(), right_atom.element()] != elements {
        return None;
    }
    component
        .bonds()
        .iter()
        .find(|candidate| {
            let [a, b] = candidate.atoms();
            (a, b) == (left, right) || (a, b) == (right, left)
        })
        .map(TemplateBond::metal_coordination)
}

fn is_hydrogen(element: Element) -> bool {
    element.atomic_number() == 1
}

fn is_carbon(element: Element) -> bool {
    element.atomic_number() == 6
}

struct Resolver<'r, 'g> {
    catalog: &'r ConnectivityCatalog,
    graph: &'g mut MoleculeEditor,
    source_aromatic_bonds: &'g mut BTreeSet<crate::molecule::BondId>,
    evidence: &'g mut BondEvidenceMap,
    issues: &'g mut Vec<raw::StructureInterpretIssue>,
    model_id: &'r str,
    cancellation: Option<&'r Cancellation>,
}

impl<'r> Resolver<'r, '_> {
    fn resolve(
        &mut self,
        residues: BTreeMap<ResidueKey, ResidueAtoms>,
        listed: &[ListedBond],
    ) -> Result<(), FormatError> {
        let mut residue_of = BTreeMap::new();
        for (step, residue) in residues.values().enumerate() {
            raw::poll_every(self.cancellation, step)?;
            for site in &residue.sites {
                residue_of.insert(site.atom, &residue.key);
            }
        }
        let mut listed_in = BTreeMap::<&ResidueKey, Vec<&ListedBond>>::new();
        let mut listed_across = Vec::new();
        for bond in listed {
            let [left, right] = bond.atoms;
            match (residue_of.get(&left), residue_of.get(&right)) {
                (Some(&left), Some(&right)) if left == right => {
                    listed_in.entry(left).or_default().push(bond);
                }
                _ => listed_across.push(bond),
            }
        }
        let mut templated = Vec::new();
        let mut stray_hydrogens = Vec::new();
        for residue in residues.values() {
            raw::poll(self.cancellation)?;
            let listed = listed_in.remove(&residue.key).unwrap_or_default();
            let component_id = normalized(&residue.key.component_id);
            let template = self.catalog.residue_template(residue);
            let force_field = template.as_ref().and_then(|template| template.force_field);
            let templated_before = templated.len();
            let hydrogens = || {
                residue
                    .sites
                    .iter()
                    .filter(|site| is_hydrogen(site.element))
                    .map(|site| (residue.key.clone(), site.atom))
            };
            if let Some(bonds) = self.catalog.component_bonds.get(&component_id) {
                // Stated orders come first, so that a table bond of another
                // order fails.
                let stated = listed
                    .iter()
                    .copied()
                    .filter(|bond| bond.order_source == OrderSource::Stated)
                    .collect::<Vec<_>>();
                self.add_listed_bonds(&stated)?;
                for bond in bonds {
                    let (Some(&left), Some(&right)) = (
                        residue.atoms.get(&bond.atom_1),
                        residue.atoms.get(&bond.atom_2),
                    ) else {
                        continue;
                    };
                    let order =
                        if self.is_coordination(bond, residue, template.as_ref(), [left, right]) {
                            SourceBondOrder::Localized(BondOrder::Dative)
                        } else {
                            bond.order
                        };
                    self.add_bond(left, right, order, raw::BondEvidence::ComponentTable)?;
                }
                self.add_listed_bonds(&listed)?;
                // Hydrogens the table leaves unbonded, such as N-terminal or
                // 5'-hydroxyl hydrogens of a polymer residue, are strays.
                stray_hydrogens.extend(hydrogens());
                // The block's own table supplies the bonds; a template whose
                // atoms the residue matches still supplies charges.
                if let Some(template) = template {
                    if let Some((matched, _)) = match_template(residue, template.component, false)?
                    {
                        templated.push(TemplatedResidue {
                            key: residue.key.clone(),
                            component: template.component,
                            variants: template.variants,
                            matched,
                            stated: Vec::new(),
                        });
                    }
                }
            } else if let Some(template) = template {
                let component = template.component;
                // A residue whose atoms fit no template is read as its
                // records bond it when they bond every heavy atom, as
                // writers state every bond of one; otherwise it fails.
                let complete = !listed.is_empty()
                    && residue
                        .sites
                        .iter()
                        .filter(|site| !is_hydrogen(site.element))
                        .all(|site| listed.iter().any(|bond| bond.atoms.contains(&site.atom)));
                if let Some((matched, strays)) = match_template(residue, component, !complete)? {
                    let stated =
                        self.add_checked_listed_bonds(residue, component, &matched, &listed)?;
                    self.add_template_bonds(component, &matched)?;
                    // Listed bonds of atoms the template does not name, such
                    // as hydrogens, follow its own.
                    self.add_listed_bonds(&listed)?;
                    stray_hydrogens
                        .extend(strays.into_iter().map(|atom| (residue.key.clone(), atom)));
                    templated.push(TemplatedResidue {
                        key: residue.key.clone(),
                        component,
                        variants: template.variants,
                        matched,
                        stated,
                    });
                } else {
                    self.issues
                        .push(raw::StructureInterpretIssue::ResidueFitsNoTemplate {
                            residue: residue.key.id(self.model_id),
                            component_id: residue.key.component_id.clone(),
                        });
                    self.add_untemplated_listed_bonds(residue, &listed)?;
                    stray_hydrogens.extend(hydrogens());
                }
            } else if !listed.is_empty() {
                // The source lists this residue's bonds; inference would only
                // guess at them. Hydrogens the listing leaves out are strays.
                self.add_untemplated_listed_bonds(residue, &listed)?;
                stray_hydrogens.extend(hydrogens());
            } else if residue.sites.len() > 1 {
                self.infer_residue_bonds(residue)?;
            }
            // Reported only when the template is applied.
            if let Some(force_field) = force_field.filter(|_| templated.len() > templated_before) {
                self.issues
                    .push(raw::StructureInterpretIssue::ForceFieldTemplate {
                        residue: residue.key.id(self.model_id),
                        component_id: residue.key.component_id.clone(),
                        force_field: force_field.to_owned(),
                    });
            }
        }
        raw::poll(self.cancellation)?;
        self.add_listed_bonds(&listed_across)?;
        raw::poll(self.cancellation)?;
        self.apply_branch_links(&residues)?;
        raw::poll(self.cancellation)?;
        self.apply_polymer_links(&residues)?;
        raw::poll(self.cancellation)?;
        self.attach_stray_hydrogens(&residues, stray_hydrogens)?;
        for residue in &templated {
            raw::poll(self.cancellation)?;
            self.assign_template_charges(residue, &residues[&residue.key])?;
        }
        let templated = templated
            .iter()
            .map(|residue| &residue.key)
            .collect::<BTreeSet<_>>();
        for residue in residues.values() {
            raw::poll(self.cancellation)?;
            if !templated.contains(&residue.key) {
                self.assign_implied_charges(residue)?;
            }
        }
        raw::poll(self.cancellation)?;
        self.report_fragmented_residues(&residues);
        Ok(())
    }

    /// Charges the atoms of a residue that no template gives charges, as
    /// their bonds imply ([`implied_charges`]). A source can state bond
    /// orders but no charges, as the component tables of wwPDB entry files
    /// do. A neutral atom with more bonds than its valence allows, such as
    /// the methylated N7 of 7-methylguanosine, takes the charge that makes
    /// its valence valid, +1 before -1, and an onium centre, such as a sulfur
    /// bonded to three carbons, takes +1 rather than a hydrogen. A singly
    /// bonded chalcogen or boron neighbor of either takes -1 with its +1, as
    /// the oxygen of DMSO drawn with a single bond. An onium centre without
    /// such a neighbor stays neutral when some of the residue's bonds have
    /// no stated order: they come from atom distances
    /// ([`raw::StructureInterpretIssue::ResidueBondsInferred`]) or from
    /// records without orders
    /// ([`raw::StructureInterpretIssue::ResidueBondOrdersUnstated`]), which
    /// then lists it. Metal donors ([`is_metal_donor`]) take no hydrogens. Stated
    /// charges are kept. Alkali, alkaline-earth, and p-block metals at a
    /// valence their element does not allow are reported.
    fn assign_implied_charges(&mut self, residue: &ResidueAtoms) -> Result<(), FormatError> {
        let id = residue.key.id(self.model_id);
        let unstated_orders = self.issues.iter().position(|issue| {
            matches!(
                issue,
                raw::StructureInterpretIssue::ResidueBondsInferred { residue, .. }
                    | raw::StructureInterpretIssue::ResidueBondOrdersUnstated { residue, .. }
                    if *residue == id
            )
        });
        let may_charge = |atom| {
            residue
                .site(atom)
                .is_some_and(|site| site.stated_charge.is_none())
        };
        let mut charged_atoms = Vec::new();
        let mut uncharged = Vec::new();
        for site in &residue.sites {
            if site.stated_charge.is_some() {
                continue;
            }
            let atom = self
                .graph
                .require_atom(site.atom)
                .map_err(interpret_error)?;
            let charges = match implied_charges(self.graph, site.atom, atom, may_charge) {
                None => continue,
                Some(ImpliedCharges::Onium) if unstated_orders.is_some() => {
                    let hydrogens = onium_hydrogens(self.graph, site.atom, atom).unwrap_or(0);
                    uncharged.push((site.name.clone(), hydrogens));
                    continue;
                }
                Some(ImpliedCharges::Onium) => vec![(site.atom, 1)],
                Some(ImpliedCharges::OverBonded(charge)) => vec![(site.atom, charge)],
                Some(ImpliedCharges::Separated(partner)) => vec![(site.atom, 1), (partner, -1)],
            };
            for (atom, charge) in charges {
                self.graph
                    .require_atom_mut(atom)
                    .map_err(interpret_error)?
                    .formal_charge = charge;
                let name = &residue.site(atom).expect("charged atoms are sites").name;
                charged_atoms.push((name.clone(), charge));
            }
        }
        let component_id = &residue.key.component_id;
        let mut donors = Vec::new();
        for site in &residue.sites {
            let atom = self
                .graph
                .require_atom(site.atom)
                .map_err(interpret_error)?;
            if is_metal_donor(self.graph, site.atom, atom) {
                self.graph
                    .require_atom_mut(site.atom)
                    .map_err(interpret_error)?
                    .hydrogens = ImplicitHydrogens::Fixed(0);
                donors.push(site.name.clone());
            }
        }
        for site in &residue.sites {
            let atom = self
                .graph
                .require_atom(site.atom)
                .map_err(interpret_error)?;
            if let Some(valence) = incomplete_metal_valence(self.graph, site.atom, atom) {
                self.issues
                    .push(raw::StructureInterpretIssue::MetalValenceIncomplete {
                        residue: id.clone(),
                        component_id: component_id.clone(),
                        atom_name: site.name.clone(),
                        valence,
                    });
            }
        }
        if !charged_atoms.is_empty() {
            self.issues
                .push(raw::StructureInterpretIssue::ChargesFromBonds {
                    residue: id.clone(),
                    component_id: component_id.clone(),
                    charged_atoms,
                });
        }
        if !donors.is_empty() {
            self.issues.push(raw::StructureInterpretIssue::MetalDonors {
                residue: id.clone(),
                component_id: component_id.clone(),
                atom_names: donors,
            });
        }
        if let Some(index) = unstated_orders {
            match &mut self.issues[index] {
                raw::StructureInterpretIssue::ResidueBondsInferred {
                    neutral_centres, ..
                }
                | raw::StructureInterpretIssue::ResidueBondOrdersUnstated {
                    neutral_centres, ..
                } => *neutral_centres = uncharged,
                _ => unreachable!("the index finds an issue about unstated orders"),
            }
        }
        Ok(())
    }

    fn add_bond(
        &mut self,
        left: AtomId,
        right: AtomId,
        order: SourceBondOrder,
        source: raw::BondEvidence,
    ) -> Result<(), FormatError> {
        add_bond_if_missing(self.graph, left, right, order, self.source_aromatic_bonds)?;
        record_evidence(self.evidence, left, right, source);
        Ok(())
    }

    /// Adds the bonds records state within a templated residue before the
    /// template's bonds. A bond between two template atoms that the template
    /// lacks is added and reported. A stated order that the template does not
    /// give is added too, and the stated template bonds are returned, to be
    /// checked against the template or variant the residue takes. The
    /// template gives the other bonds their orders.
    fn add_checked_listed_bonds(
        &mut self,
        residue: &ResidueAtoms,
        component: &Component,
        matched: &BTreeMap<usize, AtomId>,
        listed: &[&ListedBond],
    ) -> Result<Vec<StatedBond>, FormatError> {
        let index_of = matched
            .iter()
            .map(|(&index, &atom)| (atom, index))
            .collect::<BTreeMap<_, _>>();
        let mut stated = Vec::new();
        for bond in listed {
            let [left, right] = bond.atoms;
            let (Some(&a), Some(&b)) = (index_of.get(&left), index_of.get(&right)) else {
                continue;
            };
            if self
                .graph
                .require_bond_between(left, right)
                .map_err(interpret_error)?
                .is_some()
            {
                continue;
            }
            let template_bond = component.bonds().iter().find(|candidate| {
                let [x, y] = candidate.atoms();
                (x, y) == (a, b) || (x, y) == (b, a)
            });
            let names = [left, right].map(|atom| {
                residue
                    .site(atom)
                    .expect("listed atoms are residue sites")
                    .name
                    .clone()
            });
            match template_bond {
                None => self
                    .issues
                    .push(raw::StructureInterpretIssue::BondNotInTemplate {
                        residue: residue.key.id(self.model_id),
                        component_id: residue.key.component_id.clone(),
                        atoms: names,
                    }),
                Some(_) if bond.order_source != OrderSource::Stated => continue,
                Some(template_bond) => {
                    let agrees = gives_order(template_bond, bond.order);
                    stated.push(StatedBond {
                        atoms: [left, right],
                        template_atoms: [a, b],
                        names,
                        order: bond.order,
                    });
                    if agrees {
                        continue;
                    }
                }
            }
            self.add_bond(
                left,
                right,
                SourceBondOrder::Localized(bond.order),
                bond.evidence.clone(),
            )?;
        }
        Ok(stated)
    }

    /// Reports each stated order of a templated residue that `candidate`,
    /// its template or a variant, does not give.
    fn report_stated_orders(&mut self, templated: &TemplatedResidue<'_>, candidate: &Component) {
        for bond in &templated.stated {
            let Some(template_bond) = candidate_bond(templated.component, candidate, bond) else {
                continue;
            };
            if !gives_order(template_bond, bond.order) {
                self.issues
                    .push(raw::StructureInterpretIssue::BondOrderDiffersFromTemplate {
                        residue: templated.key.id(self.model_id),
                        component_id: templated.key.component_id.clone(),
                        atoms: bond.names.clone(),
                        order: bond.order,
                        template_order: template_bond_order(template_bond),
                    });
            }
        }
    }

    /// Adds the listed bonds of a residue that no template or bond table
    /// covers, and reports the residue when records state no order for some
    /// of its bonds between heavy atoms.
    fn add_untemplated_listed_bonds(
        &mut self,
        residue: &ResidueAtoms,
        listed: &[&ListedBond],
    ) -> Result<(), FormatError> {
        let heavy = |atom| {
            residue
                .site(atom)
                .is_some_and(|site| !is_hydrogen(site.element))
        };
        let mut unstated = 0;
        for bond in listed {
            let [left, right] = bond.atoms;
            if self
                .graph
                .require_bond_between(left, right)
                .map_err(interpret_error)?
                .is_some()
            {
                continue;
            }
            self.add_bond(
                left,
                right,
                SourceBondOrder::Localized(bond.order),
                bond.evidence.clone(),
            )?;
            if bond.order_source == OrderSource::Unstated && heavy(left) && heavy(right) {
                unstated += 1;
            }
        }
        if unstated > 0 {
            self.issues
                .push(raw::StructureInterpretIssue::ResidueBondOrdersUnstated {
                    residue: residue.key.id(self.model_id),
                    component_id: residue.key.component_id.clone(),
                    bonds: unstated,
                    neutral_centres: Vec::new(),
                });
        }
        Ok(())
    }

    /// Adds listed bonds between atoms not yet bonded; an existing bond keeps
    /// its order.
    fn add_listed_bonds(&mut self, listed: &[&ListedBond]) -> Result<(), FormatError> {
        for bond in listed {
            let [left, right] = bond.atoms;
            if self
                .graph
                .require_bond_between(left, right)
                .map_err(interpret_error)?
                .is_some()
            {
                continue;
            }
            self.add_bond(
                left,
                right,
                SourceBondOrder::Localized(bond.order),
                bond.evidence.clone(),
            )?;
        }
        Ok(())
    }

    /// Adds the template bonds among matched atoms. Bonds the source states
    /// explicitly outrank the template.
    fn add_template_bonds(
        &mut self,
        component: &Component,
        matched: &BTreeMap<usize, AtomId>,
    ) -> Result<(), FormatError> {
        for bond in component.bonds() {
            let [a, b] = bond.atoms();
            let (Some(&left), Some(&right)) = (matched.get(&a), matched.get(&b)) else {
                continue;
            };
            if self
                .graph
                .require_bond_between(left, right)
                .map_err(interpret_error)?
                .is_some()
            {
                continue;
            }
            self.add_bond(
                left,
                right,
                SourceBondOrder::of_template(bond),
                raw::BondEvidence::Template {
                    component_id: component.id().to_owned(),
                    variant: None,
                },
            )?;
        }
        Ok(())
    }

    /// Whether a bond of the block's own component table is metal
    /// coordination, which is dative whatever its stated order. A stated
    /// `pdbx_metal_coordination_flag` decides. Otherwise the flag comes from
    /// the component's template, or one of its variants, when it lists the
    /// bond, as charges do; without one, a bond between a metal and a
    /// non-metal is coordination, as for PDB `LINK` and `CONECT` records. A
    /// source and its PDB rendering then give the same bonds.
    fn is_coordination(
        &self,
        bond: &ComponentBond,
        residue: &ResidueAtoms,
        template: Option<&ResidueTemplate<'_>>,
        atoms: [AtomId; 2],
    ) -> bool {
        if let Some(stated) = bond.coordination {
            return stated;
        }
        let elements = atoms.map(|atom| {
            residue
                .site(atom)
                .expect("component-table atoms are residue sites")
                .element
        });
        let templated = template.and_then(|template| {
            std::iter::once(template.component)
                .chain(template.variants.iter().copied())
                .find_map(|component| template_coordination(component, bond, elements))
        });
        templated.unwrap_or_else(|| raw::is_metal(elements[0]) != raw::is_metal(elements[1]))
    }

    /// Bonds every pair of a residue's atoms within covalent bonding
    /// distance, when the caller allows it.
    fn infer_residue_bonds(&mut self, residue: &ResidueAtoms) -> Result<(), FormatError> {
        if !self.catalog.infer_bonds_from_distance {
            return Err(raw::structure_error(
                DiagnosticKind::Unsupported,
                residue.sites.first().map(|site| site.line),
                format!(
                    "Component {:?} has no template or bond table. Supply a template or turn on distance inference.",
                    residue.key.component_id
                ),
            ));
        }
        let mut bonds = 0;
        for (index, left) in residue.sites.iter().enumerate() {
            // The pair scan is quadratic in a residue's sites, so each row polls.
            raw::poll(self.cancellation)?;
            for right in &residue.sites[index + 1..] {
                if within_bonding_distance((left.element, left.point), (right.element, right.point))
                    && self
                        .graph
                        .require_bond_between(left.atom, right.atom)
                        .map_err(interpret_error)?
                        .is_none()
                {
                    self.add_bond(
                        left.atom,
                        right.atom,
                        SourceBondOrder::Localized(BondOrder::Single),
                        raw::BondEvidence::Distance,
                    )?;
                    bonds += 1;
                }
            }
        }
        self.issues
            .push(raw::StructureInterpretIssue::ResidueBondsInferred {
                residue: residue.key.id(self.model_id),
                component_id: residue.key.component_id.clone(),
                bonds,
                neutral_centres: Vec::new(),
            });
        Ok(())
    }

    /// Links the residues of each copy of a branched entity as its
    /// `_pdbx_entity_branch_link` rows state: both ends lie in the same asym,
    /// which `_pdbx_branch_scheme` lists, and the linking atoms within bonding
    /// distance. A longer gap is a reported break, as between polymer
    /// residues.
    fn apply_branch_links(
        &mut self,
        residues: &BTreeMap<ResidueKey, ResidueAtoms>,
    ) -> Result<(), FormatError> {
        for link in &self.catalog.branch_links {
            raw::poll(self.cancellation)?;
            let asyms = self
                .catalog
                .branch_sites
                .iter()
                .filter(|site| site.entity_id == link.entity_id)
                .map(|site| site.asym_id.as_str())
                .collect::<BTreeSet<_>>();
            for asym_id in asyms {
                let resolve = |number| {
                    resolve_branch_residue(residues, self.catalog, &link.entity_id, asym_id, number)
                };
                let (Some(left), Some(right)) = (resolve(link.number_1), resolve(link.number_2))
                else {
                    continue;
                };
                let (Some(&left_atom), Some(&right_atom)) =
                    (left.atoms.get(&link.atom_1), right.atoms.get(&link.atom_2))
                else {
                    continue;
                };
                let (left_site, right_site) = (
                    left.site(left_atom)
                        .expect("linking atom is a residue site"),
                    right
                        .site(right_atom)
                        .expect("linking atom is a residue site"),
                );
                if within_bonding_distance(
                    (left_site.element, left_site.point),
                    (right_site.element, right_site.point),
                ) {
                    self.add_bond(
                        left_atom,
                        right_atom,
                        link.order,
                        raw::BondEvidence::BranchLink,
                    )?;
                } else {
                    self.issues.push(raw::StructureInterpretIssue::ChainBreak {
                        residue: left.key.id(self.model_id),
                        next: right.key.id(self.model_id),
                    });
                }
            }
        }
        Ok(())
    }

    /// Links residues adjacent in sequence within one chain: peptide `C`-`N`
    /// or nucleotide `O3'`-`P`, accepted only within bonding distance. The
    /// polymer type comes from `_entity_poly`, or from the residues'
    /// templates when the source does not describe the entity. A link one of
    /// whose atoms a connection record bonds to an atom of another residue,
    /// as in an isopeptide, is not made, and the report names it. Neither is
    /// one whose atoms are missing; the report names it when each residue
    /// has a linking atom and no bond joins them, so that models without
    /// backbone atoms, such as CA-only or coarse-grained ones, and runs of
    /// waters and ions are not reported. Entities whose `_entity_poly`
    /// states a nonstandard linkage take the same rule: the flag does not
    /// say which links differ.
    fn apply_polymer_links(
        &mut self,
        residues: &BTreeMap<ResidueKey, ResidueAtoms>,
    ) -> Result<(), FormatError> {
        let mut chains = BTreeMap::<(String, Option<String>), Vec<&ResidueAtoms>>::new();
        for residue in residues.values() {
            if residue.key.polymer_position.is_none() {
                continue;
            }
            chains
                .entry((residue.key.asym_id.clone(), residue.key.entity_id.clone()))
                .or_default()
                .push(residue);
        }
        for ((_, entity), chain) in &mut chains {
            chain.sort_by_key(|residue| residue.key.polymer_position);
            let declared = entity
                .as_ref()
                .and_then(|entity| self.catalog.polymer_types.get(entity));
            for (step, pair) in chain.windows(2).enumerate() {
                raw::poll_every(self.cancellation, step)?;
                let (left, right) = (pair[0], pair[1]);
                let (Some(left_seq), Some(right_seq)) =
                    (left.key.polymer_position, right.key.polymer_position)
                else {
                    continue;
                };
                // Positions count unmodelled residues, as mmCIF label_seq_id
                // does, so only adjacent values prove that the two observed
                // residues are consecutive.
                if right_seq != left_seq.saturating_add(1) {
                    continue;
                }
                let link = match declared {
                    Some(kind) if is_polypeptide(kind) => Some(PolymerLink::Peptide),
                    Some(kind) if is_nucleic_acid(kind) => Some(PolymerLink::Nucleotide),
                    Some(_) => None,
                    None => PolymerLink::undeclared(
                        self.catalog.linking(&left.key.component_id),
                        self.catalog.linking(&right.key.component_id),
                    ),
                };
                let Some(link) = link else {
                    continue;
                };
                let site = |atom| {
                    left.site(atom)
                        .or_else(|| right.site(atom))
                        .expect("linking atom is a residue site")
                };
                let outcome = standard_link(
                    link,
                    |name| left.atoms.get(name).copied(),
                    |name| right.atoms.get(name).copied(),
                    |atom, partner| {
                        let residue = if left.site(atom).is_some() {
                            left
                        } else {
                            right
                        };
                        self.bonded_elsewhere(atom, residue, partner)
                    },
                    |left_atom, right_atom| {
                        let (left_site, right_site) = (site(left_atom), site(right_atom));
                        within_bonding_distance(
                            (left_site.element, left_site.point),
                            (right_site.element, right_site.point),
                        )
                    },
                )?;
                let issue = match outcome {
                    StandardLink::Linked([left_atom, right_atom]) => {
                        self.add_bond(
                            left_atom,
                            right_atom,
                            SourceBondOrder::Localized(BondOrder::Single),
                            raw::BondEvidence::PolymerLink,
                        )?;
                        continue;
                    }
                    StandardLink::Missing => {
                        if !has_linking_atom(left)
                            || !has_linking_atom(right)
                            || self.residues_bonded(left, right)?
                        {
                            continue;
                        }
                        raw::StructureInterpretIssue::LinkingAtomMissing {
                            residue: left.key.id(self.model_id),
                            next: right.key.id(self.model_id),
                        }
                    }
                    StandardLink::Skipped => raw::StructureInterpretIssue::PolymerLinkSkipped {
                        residue: left.key.id(self.model_id),
                        next: right.key.id(self.model_id),
                    },
                    StandardLink::TooFar => raw::StructureInterpretIssue::ChainBreak {
                        residue: left.key.id(self.model_id),
                        next: right.key.id(self.model_id),
                    },
                };
                self.issues.push(issue);
            }
        }
        Ok(())
    }

    /// Whether a bond joins an atom of `left` to one of `right`.
    fn residues_bonded(
        &self,
        left: &ResidueAtoms,
        right: &ResidueAtoms,
    ) -> Result<bool, FormatError> {
        for site in &left.sites {
            if self
                .graph
                .require_neighbors(site.atom)
                .map_err(interpret_error)?
                .any(|neighbor| right.site(neighbor).is_some())
            {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Whether a bond joins `atom` of `residue` to an atom of another residue
    /// than `partner`.
    fn bonded_elsewhere(
        &self,
        atom: AtomId,
        residue: &ResidueAtoms,
        partner: AtomId,
    ) -> Result<bool, FormatError> {
        Ok(self
            .graph
            .require_neighbors(atom)
            .map_err(interpret_error)?
            .any(|neighbor| neighbor != partner && residue.site(neighbor).is_none()))
    }

    /// Bonds each unbonded hydrogen that its template or listed bonds do not
    /// cover to the nearest heavy atom of its residue within bonding distance.
    fn attach_stray_hydrogens(
        &mut self,
        residues: &BTreeMap<ResidueKey, ResidueAtoms>,
        stray_hydrogens: Vec<(ResidueKey, AtomId)>,
    ) -> Result<(), FormatError> {
        for (key, hydrogen) in stray_hydrogens {
            // Each hydrogen scans its residue's sites.
            raw::poll(self.cancellation)?;
            let residue = &residues[&key];
            let site = residue
                .site(hydrogen)
                .expect("stray hydrogen is a residue site");
            if self
                .graph
                .require_incident_bonds(hydrogen)
                .map_err(interpret_error)?
                .next()
                .is_some()
            {
                continue;
            }
            let distance = |other: &ResidueSite| {
                (other.point.x - site.point.x).powi(2)
                    + (other.point.y - site.point.y).powi(2)
                    + (other.point.z - site.point.z).powi(2)
            };
            let nearest = residue
                .sites
                .iter()
                .filter(|other| !is_hydrogen(other.element))
                .filter(|other| {
                    within_bonding_distance(
                        (site.element, site.point),
                        (other.element, other.point),
                    )
                })
                .min_by(|a, b| distance(a).total_cmp(&distance(b)));
            let Some(heavy) = nearest else {
                continue;
            };
            self.add_bond(
                hydrogen,
                heavy.atom,
                SourceBondOrder::Localized(BondOrder::Single),
                raw::BondEvidence::HydrogenDistance,
            )?;
            self.issues
                .push(raw::StructureInterpretIssue::HydrogenAttachedByDistance {
                    residue: key.id(self.model_id),
                    atom_name: site.name.clone(),
                    heavy_atom_name: heavy.name.clone(),
                });
        }
        Ok(())
    }

    /// Sets formal charges of a templated residue. Stated charges win, and
    /// only the template or a variant that gives every stated bond order is
    /// a candidate; without one, charges follow from valence, as for a
    /// residue without a template. Without model hydrogens, the template's
    /// charges apply. With them, the template
    /// or protonation variant matching the hydrogens on every heteroatom is
    /// preferred; otherwise one whose differences lie on heteroatoms with only
    /// single bonds, each added or removed hydrogen changing the charge by one.
    /// Hydrogens absent from carbon are left to valence.
    fn assign_template_charges(
        &mut self,
        templated: &TemplatedResidue<'_>,
        residue: &ResidueAtoms,
    ) -> Result<(), FormatError> {
        let named = templated
            .matched
            .iter()
            .map(|(&index, &atom)| (templated.component.atoms()[index].name(), atom))
            .collect::<Vec<_>>();
        let has_hydrogens = residue.sites.iter().any(|site| is_hydrogen(site.element));
        // Stated orders outrank hydrogens: only a template or variant that
        // gives every stated order is a candidate.
        let exact = |charges: &[TemplateCharge]| charges.iter().all(|charge| charge.delta == 0);
        let mut evaluated = Vec::new();
        let mut excluded_exact = false;
        for candidate in
            std::iter::once(templated.component).chain(templated.variants.iter().copied())
        {
            let Some(charges) =
                self.candidate_charges(candidate, &named, residue, has_hydrogens)?
            else {
                continue;
            };
            if templated
                .stated
                .iter()
                .all(|bond| candidate_gives(templated.component, candidate, bond))
            {
                evaluated.push((candidate, charges));
            } else {
                excluded_exact |= exact(&charges);
            }
        }
        if evaluated.is_empty() && !templated.stated.is_empty() {
            // No candidate fits both the stated orders and the hydrogens:
            // the stated orders stay, and charges follow from valence. A
            // residue without hydrogens has nothing they disagree with.
            self.report_stated_orders(templated, templated.component);
            if has_hydrogens {
                self.report_disagreement(templated);
            }
            return self.assign_implied_charges(residue);
        }
        let chosen = evaluated
            .iter()
            .find(|(_, charges)| exact(charges))
            .or_else(|| evaluated.first())
            .ok_or_else(|| {
                raw::structure_error(
                    DiagnosticKind::Inconsistent,
                    residue.sites.first().map(|site| site.line),
                    format!(
                        "The hydrogens of a {:?} residue match neither its template nor a protonation variant.",
                        templated.component.id()
                    ),
                )
            })?;
        let (candidate, charges) = chosen;
        let mut charged_atoms = Vec::new();
        for charge in charges {
            let site = residue
                .site(charge.atom)
                .expect("charged atom is a residue site");
            let value = site.stated_charge.unwrap_or(charge.value);
            self.graph
                .require_atom_mut(charge.atom)
                .map_err(interpret_error)?
                .formal_charge = value;
            if charge.delta != 0 {
                charged_atoms.push((site.name.clone(), value));
            }
        }
        self.report_stated_orders(templated, candidate);
        if candidate.variant().is_some() {
            let stated = templated
                .stated
                .iter()
                .map(|bond| bond.atoms)
                .collect::<Vec<_>>();
            self.apply_variant_bonds(candidate, &named, &stated)?;
        }
        // A candidate the stated orders exclude fits the hydrogens without
        // charges, and the chosen one does not.
        let disagree = excluded_exact && !exact(charges);
        if candidate.variant().is_some() || !charged_atoms.is_empty() {
            self.issues
                .push(raw::StructureInterpretIssue::ProtonationFromHydrogens {
                    residue: templated.key.id(self.model_id),
                    component_id: templated.component.id().to_owned(),
                    variant: candidate.variant().map(str::to_owned),
                    charged_atoms,
                });
        }
        if disagree {
            self.report_disagreement(templated);
        }
        Ok(())
    }

    /// Reports a residue whose stated bond orders no candidate gives with
    /// its hydrogens.
    fn report_disagreement(&mut self, templated: &TemplatedResidue<'_>) {
        self.issues.push(
            raw::StructureInterpretIssue::StatedOrdersDisagreeWithHydrogens {
                residue: templated.key.id(self.model_id),
                component_id: templated.component.id().to_owned(),
            },
        );
    }

    /// The charges a candidate template gives the residue's atoms, or `None`
    /// when the model's hydrogens contradict it.
    fn candidate_charges(
        &self,
        candidate: &Component,
        named: &[(&str, AtomId)],
        residue: &ResidueAtoms,
        has_hydrogens: bool,
    ) -> Result<Option<Vec<TemplateCharge>>, FormatError> {
        let differences = hydrogen_differences(
            candidate,
            named,
            has_hydrogens,
            |atom| is_hydrogen(residue.site(atom).expect("matched site").element),
            |atom| {
                let present = self
                    .graph
                    .require_neighbors(atom)
                    .map_err(interpret_error)?
                    .filter(|&neighbor| {
                        self.graph
                            .atom(neighbor)
                            .is_some_and(|neighbor| is_hydrogen(neighbor.element))
                    })
                    .count();
                Ok((present, self.links_outside(atom, residue)?))
            },
        )?;
        let Some(differences) = differences else {
            return Ok(None);
        };
        differences
            .into_iter()
            .map(|(atom, index, difference)| {
                let base = candidate.atoms()[index].formal_charge().unwrap_or(0);
                let delta = i8::try_from(difference).map_err(|_| {
                    chemistry_error("A residue's hydrogen count differs from its template by more than a formal charge holds.")
                })?;
                Ok(TemplateCharge {
                    atom,
                    value: base.saturating_add(delta),
                    delta,
                })
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Some)
    }

    /// The valence that bonds to atoms of other residues give a residue
    /// atom (see [`link_valence`]).
    fn links_outside(&self, atom: AtomId, residue: &ResidueAtoms) -> Result<usize, FormatError> {
        Ok(self
            .graph
            .require_incident_bonds(atom)
            .map_err(interpret_error)?
            .filter(|(_, bond)| residue.site(bond.other_atom(atom)).is_none())
            .map(|(_, bond)| link_valence(bond.order()))
            .sum())
    }

    /// Replaces template bond orders with a variant's, except the orders
    /// records state between the `stated` atom pairs.
    fn apply_variant_bonds(
        &mut self,
        variant: &Component,
        named: &[(&str, AtomId)],
        stated: &[[AtomId; 2]],
    ) -> Result<(), FormatError> {
        let atom_of = |index: usize| {
            let name = variant.atoms()[index].name();
            named
                .iter()
                .find(|(candidate, _)| *candidate == name)
                .map(|&(_, atom)| atom)
        };
        for bond in variant.bonds() {
            let [a, b] = bond.atoms();
            let (Some(left), Some(right)) = (atom_of(a), atom_of(b)) else {
                continue;
            };
            if stated.contains(&[left, right]) || stated.contains(&[right, left]) {
                continue;
            }
            let Some(existing) = self
                .graph
                .require_bond_between(left, right)
                .map_err(interpret_error)?
            else {
                continue;
            };
            self.graph
                .require_bond_mut(existing)
                .map_err(interpret_error)?
                .set_order(SourceBondOrder::of_template(bond).staged_order());
            self.evidence.insert(
                (left.min(right), left.max(right)),
                raw::BondEvidence::Template {
                    component_id: variant.id().to_owned(),
                    variant: variant.variant().map(str::to_owned),
                },
            );
        }
        Ok(())
    }

    fn report_fragmented_residues(&mut self, residues: &BTreeMap<ResidueKey, ResidueAtoms>) {
        let components = self.graph.connected_components();
        let mut component_of = vec![None; self.graph.atom_slot_count()];
        for (index, component) in components.iter().enumerate() {
            for atom in component {
                component_of[atom.index()] = Some(index);
            }
        }
        for residue in residues.values() {
            let fragments = residue
                .sites
                .iter()
                .filter_map(|site| component_of.get(site.atom.index()).copied().flatten())
                .collect::<BTreeSet<_>>()
                .len();
            if fragments > 1 {
                self.issues
                    .push(raw::StructureInterpretIssue::ResidueFragmented {
                        residue: residue.key.id(self.model_id),
                        fragments,
                    });
            }
        }
    }
}

/// Template atom indices matched to staged atoms, and unnamed hydrogens.
type TemplateMatch = (BTreeMap<usize, AtomId>, Vec<AtomId>);

/// Matches residue atoms to template atoms by name and element, returning the
/// match and the hydrogens the template does not name. A strict match fails on
/// any other unknown or mismatched atom; a lenient one returns `None`.
fn match_template(
    residue: &ResidueAtoms,
    component: &Component,
    strict: bool,
) -> Result<Option<TemplateMatch>, FormatError> {
    let mut matched = BTreeMap::new();
    let mut strays = Vec::new();
    for site in &residue.sites {
        let failure = |message: String| {
            if strict {
                Err(raw::structure_error(
                    DiagnosticKind::Inconsistent,
                    Some(site.line),
                    message,
                ))
            } else {
                Ok(None)
            }
        };
        match component.atom(&site.name) {
            Some((index, template_atom)) => {
                if template_atom.element() != site.element {
                    return failure(format!(
                        "Atom {:?} of component {:?} is {} in the model but {} in its template.",
                        site.name,
                        component.id(),
                        site.element.symbol(),
                        template_atom.element().symbol()
                    ));
                }
                if matched.insert(index, site.atom).is_some() {
                    return failure(format!(
                        "Two atoms of one {:?} residue match template atom {:?}.",
                        component.id(),
                        template_atom.name()
                    ));
                }
            }
            None if is_hydrogen(site.element) => strays.push(site.atom),
            None => {
                return failure(format!(
                    "Atom {:?} is not in component {:?}.",
                    site.name,
                    component.id()
                ))
            }
        }
    }
    Ok(Some((matched, strays)))
}

/// One residue atom, its template atom index, and how many hydrogens it has
/// beyond the template's.
pub(crate) type HydrogenDifference<A> = (A, usize, isize);

/// How many hydrogens each named atom of a residue has beyond those
/// `candidate` gives it: `(atom, template index, difference)`, or `None`
/// when the model's hydrogens contradict the candidate. `named` pairs
/// template atom names with residue atoms. Differences are counted only on
/// heteroatoms of a residue with model hydrogens, and only where they can be
/// charges: on an atom with only single template bonds, or by one on a ring
/// nitrogen. `hydrogens` gives an atom's bonded hydrogens and the valence
/// its bonds to other residues give it, which takes the place of hydrogens.
pub(crate) fn hydrogen_differences<A: Copy, E>(
    candidate: &Component,
    named: &[(&str, A)],
    has_hydrogens: bool,
    is_hydrogen_atom: impl Fn(A) -> bool,
    mut hydrogens: impl FnMut(A) -> Result<(usize, usize), E>,
) -> Result<Option<Vec<HydrogenDifference<A>>>, E> {
    let mut differences = Vec::new();
    for &(name, atom) in named {
        let Some((index, template_atom)) = candidate.atom(name) else {
            // Variants only remove hydrogens; a model hydrogen the variant
            // lacks shows up on its heavy atom below.
            if is_hydrogen_atom(atom) {
                continue;
            }
            return Ok(None);
        };
        let element = template_atom.element();
        let mut difference = 0;
        if has_hydrogens && !is_hydrogen(element) && !is_carbon(element) {
            let (present, links) = hydrogens(atom)?;
            let expected = expected_hydrogens(candidate, index, links);
            if present != expected {
                let only_single = candidate
                    .bonds()
                    .iter()
                    .filter(|bond| bond.atoms().contains(&index))
                    .all(|bond| {
                        bond.order() == TemplateOrder::Single
                            && !bond.aromatic()
                            && !bond.metal_coordination()
                    });
                difference = present as isize - expected as isize;
                // A ring nitrogen takes or loses one proton whatever its
                // bond orders, as a protonated cytosine N3 or a
                // deprotonated imidazole nitrogen does.
                let ring_nitrogen = element.atomic_number() == 7
                    && difference.abs() == 1
                    && in_ring(candidate, index);
                if !only_single && !ring_nitrogen {
                    return Ok(None);
                }
            }
        }
        differences.push((atom, index, difference));
    }
    Ok(Some(differences))
}

/// A formal charge from a template, after hydrogen differences.
struct TemplateCharge {
    atom: AtomId,
    value: i8,
    delta: i8,
}

/// Whether a template atom lies on a ring: some neighbour reaches it again
/// without the bond between them.
fn in_ring(component: &Component, index: usize) -> bool {
    let neighbours = |atom: usize| {
        component.bonds().iter().filter_map(move |bond| {
            let [a, b] = bond.atoms();
            match (a == atom, b == atom) {
                (true, _) => Some(b),
                (_, true) => Some(a),
                _ => None,
            }
        })
    };
    neighbours(index).any(|start| {
        let mut seen = BTreeSet::from([index, start]);
        let mut stack = vec![start];
        while let Some(atom) = stack.pop() {
            for next in neighbours(atom) {
                if next == index && atom != start {
                    return true;
                }
                if seen.insert(next) {
                    stack.push(next);
                }
            }
        }
        false
    })
}

/// The valence a bond to another residue gives an atom, which is the number
/// of hydrogens it replaces: one for a single bond, two for a double bond.
/// Zero-order and dative bonds replace none.
pub(crate) const fn link_valence(order: BondOrder) -> usize {
    match order {
        BondOrder::Zero | BondOrder::Dative => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
    }
}

/// The hydrogens a template atom carries when bonds to other residues give
/// it `links` valence. Each unit takes the place of one hydrogen, and an
/// atom that links loses all its leaving hydrogens: a peptide bond replaces
/// the amine's leaving `H2`, a disulfide replaces cysteine's `HG`, which the
/// dictionary does not mark leaving, and a Schiff base's double bond
/// replaces two of lysine's `NZ` hydrogens.
fn expected_hydrogens(component: &Component, index: usize, links: usize) -> usize {
    let (mut hydrogens, mut leaving) = (0, 0);
    for bond in component.bonds() {
        let [a, b] = bond.atoms();
        let other = match (a == index, b == index) {
            (true, _) => b,
            (_, true) => a,
            _ => continue,
        };
        let atom = &component.atoms()[other];
        if is_hydrogen(atom.element()) {
            hydrogens += 1;
            leaving += usize::from(atom.leaving());
        }
    }
    if links == 0 {
        hydrogens
    } else {
        hydrogens.saturating_sub(leaving.max(links))
    }
}

/// The residue of one asym that a branch scheme row numbers, if modelled.
fn resolve_branch_residue<'a>(
    residues: &'a BTreeMap<ResidueKey, ResidueAtoms>,
    catalog: &ConnectivityCatalog,
    entity_id: &str,
    asym_id: &str,
    number: i32,
) -> Option<&'a ResidueAtoms> {
    catalog
        .branch_sites
        .iter()
        .filter(|site| {
            site.entity_id == entity_id && site.asym_id == asym_id && site.number == number
        })
        .find_map(|site| {
            residues.values().find(|residue| {
                residue.key.asym_id == site.asym_id
                    && residue.key.entity_id.as_deref() == Some(site.entity_id.as_str())
                    && normalized(&residue.key.component_id) == normalized(&site.component_id)
                    && site.author_sequence_id.as_deref().is_none_or(|expected| {
                        residue.key.author_sequence_id.as_deref() == Some(expected)
                    })
            })
        })
}

fn add_bond_if_missing(
    graph: &mut MoleculeEditor,
    left: AtomId,
    right: AtomId,
    source_order: SourceBondOrder,
    source_aromatic_bonds: &mut BTreeSet<crate::molecule::BondId>,
) -> Result<(), FormatError> {
    let existing = graph
        .require_bond_between(left, right)
        .map_err(interpret_error)?;
    let Some(existing) = existing else {
        let bond = graph
            .add_bond(left, right, source_order.staged_order())
            .map_err(interpret_error)?;
        if source_order == SourceBondOrder::Aromatic {
            source_aromatic_bonds.insert(bond);
        }
        return Ok(());
    };
    let existing_order = graph.require_bond(existing).map_err(interpret_error)?.order;
    let consistent = match source_order {
        SourceBondOrder::Aromatic => source_aromatic_bonds.contains(&existing),
        SourceBondOrder::Localized(order) => {
            !source_aromatic_bonds.contains(&existing) && existing_order == order
        }
    };
    if !consistent {
        let source_order = match source_order {
            SourceBondOrder::Localized(order) => order.to_string(),
            SourceBondOrder::Aromatic => "aromatic".to_owned(),
        };
        return Err(chemistry_error(format!(
            "The bond records give one atom pair two orders: {existing_order} and {source_order}."
        )));
    }
    Ok(())
}

fn component_bond_order(
    value: &str,
    table: &MmcifCategory<'_>,
    row: usize,
) -> Result<SourceBondOrder, FormatError> {
    match value.to_ascii_lowercase().as_str() {
        "sing" | "poly" => Ok(SourceBondOrder::Localized(BondOrder::Single)),
        "doub" | "pi" => Ok(SourceBondOrder::Localized(BondOrder::Double)),
        "trip" => Ok(SourceBondOrder::Localized(BondOrder::Triple)),
        "quad" => Ok(SourceBondOrder::Localized(BondOrder::Quadruple)),
        "arom" | "delo" => Ok(SourceBondOrder::Aromatic),
        other => Err(row_error(
            table,
            row,
            DiagnosticKind::Unsupported,
            format!("The component bond order {other:?} is not supported."),
        )),
    }
}

fn is_yes(value: &str) -> bool {
    matches!(value.to_ascii_lowercase().as_str(), "y" | "yes")
}

fn is_polypeptide(value: &str) -> bool {
    value.to_ascii_lowercase().contains("polypeptide")
}

fn is_nucleic_acid(value: &str) -> bool {
    let value = value.to_ascii_lowercase();
    value.contains("polyribonucleotide") || value.contains("polydeoxyribonucleotide")
}

fn normalized(value: &str) -> String {
    value.to_ascii_uppercase()
}

fn required<'a>(
    table: &'a MmcifCategory<'_>,
    row: usize,
    tag: &str,
) -> Result<&'a str, FormatError> {
    let value = table.value(row, tag).ok_or_else(|| {
        row_error(
            table,
            row,
            DiagnosticKind::InvalidValue,
            format!("The required {tag} is missing."),
        )
    })?;
    value.optional_text().ok_or_else(|| {
        row_error(
            table,
            row,
            DiagnosticKind::InvalidValue,
            format!("The required {tag} is missing."),
        )
    })
}

fn required_i32(table: &MmcifCategory<'_>, row: usize, tag: &str) -> Result<i32, FormatError> {
    required(table, row, tag)?.parse::<i32>().map_err(|_| {
        row_error(
            table,
            row,
            DiagnosticKind::InvalidValue,
            format!("The {tag} value is not a whole number."),
        )
    })
}

fn optional<'a>(table: &'a MmcifCategory<'_>, row: usize, tag: &str) -> Option<&'a str> {
    table
        .value(row, tag)
        .and_then(|value| value.optional_text())
}

fn row_error(
    table: &MmcifCategory<'_>,
    row: usize,
    kind: DiagnosticKind,
    message: impl Into<String>,
) -> FormatError {
    raw::structure_error(kind, table.row_line(row), message)
}

/// Molecule construction failed with `source` as the cause.
fn interpret_error(source: impl std::error::Error + Send + Sync + 'static) -> FormatError {
    chemistry_error("Cannot build the molecule.").with_source(source)
}

fn chemistry_error(message: impl Into<String>) -> FormatError {
    raw::structure_error(DiagnosticKind::Chemistry, None, message)
}

#[cfg(test)]
mod tests;
