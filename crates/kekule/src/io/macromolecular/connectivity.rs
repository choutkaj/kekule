//! Resolves the bonds a macromolecular source does not state, from ranked
//! evidence: the source's own component bond table or connection records,
//! component templates, polymer and branch linkage, and, for residues nothing
//! else covers, covalent radii.
use crate::diagnostics::{Cancellation, DiagnosticKind, FormatError};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use crate::algorithms::{rdkit_inferred_hydrogen_count, rdkit_valence_violation, ValenceViolation};
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

mod protonation;
pub(crate) use protonation::{choose as choose_protonation, AtomEvidence, ResidueEvidence};

/// The atoms of a nucleotide's 5' phosphate; a residue without them all is a
/// 5'-hydroxyl terminus.
const FIVE_PRIME_PHOSPHATE: [&str; 6] = ["P", "OP1", "OP2", "OP3", "O1P", "O2P"];

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

/// A bond a connectivity-only record, such as PDB `CONECT`, lists between two
/// atoms of one residue. It supplements the residue's template, whose orders
/// win, and replaces distance inference for residues without one.
#[derive(Debug, Clone)]
pub(crate) struct ListedBond {
    pub(crate) atoms: [AtomId; 2],
    pub(crate) order: BondOrder,
    pub(crate) evidence: raw::BondEvidence,
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
    coordination: &BTreeMap<AtomId, usize>,
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
        coordination,
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
    /// The block's `_chem_comp_atom` charges by component and atom name;
    /// `None` where the block states the charge as unknown.
    component_charges: BTreeMap<String, BTreeMap<String, Option<i8>>>,
    polymer_types: BTreeMap<String, String>,
    nonstandard_polymer_linkage: BTreeSet<String>,
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
        if let Some(table) = block.category("_chem_comp_atom")? {
            let table = &table;
            for row in 0..table.row_count() {
                let comp_id = required(table, row, "_chem_comp_atom.comp_id")?;
                let atom = required(table, row, "_chem_comp_atom.atom_id")?;
                let Some(value) = table.value(row, "_chem_comp_atom.charge") else {
                    continue;
                };
                let charge = match value.optional_text() {
                    None => None,
                    Some(text) => Some(text.parse::<i8>().map_err(|_| {
                        row_error(
                            table,
                            row,
                            DiagnosticKind::InvalidValue,
                            format!("The component atom charge {text:?} is not a whole number."),
                        )
                    })?),
                };
                catalog
                    .component_charges
                    .entry(normalized(comp_id))
                    .or_default()
                    .insert(atom.to_owned(), charge);
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
                if optional(table, row, "_entity_poly.nstd_linkage").is_some_and(is_yes) {
                    catalog
                        .nonstandard_polymer_linkage
                        .insert(entity.to_owned());
                }
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

/// A residue whose bonds came from a template; `matched` maps template atom
/// indices to staged atoms.
struct TemplatedResidue<'c> {
    key: ResidueKey,
    component: &'c Component,
    variants: Vec<&'c Component>,
    matched: BTreeMap<usize, AtomId>,
    /// Whether the residue takes the usual pH 7 state without hydrogens: an
    /// amino acid or nucleotide read with a dictionary template.
    usual_state: bool,
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

/// The atoms a terminal hydrogen's conventional name places it on, in order
/// of preference, for hydrogens no template names: `H1`-`H3` (PDB),
/// `HT1`-`HT3` (CHARMM), and `HN1`-`HN3` on a peptide's `N`; `HO5'` and
/// `H5T` (AMBER) on a nucleotide's `O5'`, and `H3T` on its `O3'`.
fn terminal_hydrogen_anchors(linking: ComponentLinking, name: &str) -> &'static [&'static str] {
    match linking {
        ComponentLinking::Peptide
            if matches!(
                name,
                "H1" | "H2" | "H3" | "HT1" | "HT2" | "HT3" | "HN1" | "HN2" | "HN3"
            ) =>
        {
            &["N"]
        }
        ComponentLinking::Nucleotide if matches!(name, "HO5'" | "HO5*" | "H5T") => &["O5'", "O5*"],
        ComponentLinking::Nucleotide if matches!(name, "HO3*" | "H3T") => &["O3'", "O3*"],
        _ => &[],
    }
}

/// The implicit hydrogens readers fix on a residue atom of `component`
/// named `atom_name`. An atom that anchors a leaving heavy atom the residue
/// lacks (`has_atom` says which template atoms it has), and that bonds to
/// no other residue (`linked`), takes no hydrogen in that atom's place: it
/// keeps its template's other hydrogens, less its `explicit` ones. A
/// carboxyl terminus without `OXT` stays an open carbonyl, not an aldehyde.
/// `None` for an atom whose count readers infer.
pub(crate) fn fixed_hydrogens(
    component: &Component,
    atom_name: &str,
    has_atom: impl Fn(&crate::components::ComponentAtom) -> bool,
    linked: bool,
    explicit: usize,
) -> Option<u8> {
    let (index, _) = component.atom(atom_name)?;
    let anchors_missing = component.bonds().iter().any(|bond| {
        let [a, b] = bond.atoms();
        let other = match (a == index, b == index) {
            (true, _) => b,
            (_, true) => a,
            _ => return false,
        };
        let atom = &component.atoms()[other];
        atom.leaving() && !is_hydrogen(atom.element()) && !has_atom(atom)
    });
    (anchors_missing && !linked).then(|| {
        let hydrogens = protonation::expected_hydrogens(component, index, 0);
        u8::try_from(hydrogens.saturating_sub(explicit)).unwrap_or(u8::MAX)
    })
}

/// Whether a residue takes the usual pH 7 state when the model has no
/// hydrogens to decide by: an amino acid or nucleotide read with a
/// dictionary template. Force-field residue names state their protonation.
pub(crate) fn usual_state(component: &Component) -> bool {
    component.force_field().is_none()
        && matches!(
            component.linking(),
            ComponentLinking::Peptide | ComponentLinking::Nucleotide
        )
}

fn is_hydrogen(element: Element) -> bool {
    element.atomic_number() == 1
}

struct Resolver<'r, 'g> {
    catalog: &'r ConnectivityCatalog,
    graph: &'g mut MoleculeEditor,
    source_aromatic_bonds: &'g mut BTreeSet<crate::molecule::BondId>,
    evidence: &'g mut BondEvidenceMap,
    issues: &'g mut Vec<raw::StructureInterpretIssue>,
    /// Metal coordination records between residues, counted per atom.
    coordination: &'r BTreeMap<AtomId, usize>,
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
            let mut force_field = None;
            let template = self.catalog.template(&component_id).or_else(|| {
                let atoms = residue
                    .sites
                    .iter()
                    .map(|site| (site.name.as_str(), site.element));
                let (component, name) =
                    ComponentLibrary::force_field_template(&component_id, atoms)?;
                force_field = Some(name);
                Some((component, Vec::new()))
            });
            let templated_before = templated.len();
            if let Some(bonds) = self.catalog.component_bonds.get(&component_id) {
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
                stray_hydrogens.extend(
                    residue
                        .sites
                        .iter()
                        .filter(|site| is_hydrogen(site.element))
                        .map(|site| (residue.key.clone(), site.atom)),
                );
                // The block's own table supplies the bonds; a template whose
                // atoms the residue matches still supplies charges.
                if let Some((component, variants)) = template {
                    if let Some((matched, _)) = match_template(residue, component, false)? {
                        templated.push(TemplatedResidue {
                            key: residue.key.clone(),
                            component,
                            variants,
                            matched,
                            usual_state: usual_state(component),
                        });
                    }
                }
            } else if self.has_stated_internal_bond(residue)? {
                // The source enumerates this residue's bonds, as files that
                // write complete connectivity do; templates and inference
                // only fill residues the source leaves without bonds.
                self.add_listed_bonds(&listed)?;
            } else if let Some((component, variants)) = template {
                let (matched, strays) =
                    match_template(residue, component, true)?.expect("strict matches fail loudly");
                self.add_template_bonds(component, &matched)?;
                // Listed bonds the template lacks, such as those of hydrogens
                // it does not name, follow its own.
                self.add_listed_bonds(&listed)?;
                stray_hydrogens.extend(strays.into_iter().map(|atom| (residue.key.clone(), atom)));
                templated.push(TemplatedResidue {
                    key: residue.key.clone(),
                    component,
                    variants,
                    matched,
                    usual_state: usual_state(component),
                });
            } else if !listed.is_empty() {
                // The source lists this residue's bonds; inference would only
                // guess at them. Hydrogens the listing leaves out are strays.
                self.add_listed_bonds(&listed)?;
                stray_hydrogens.extend(
                    residue
                        .sites
                        .iter()
                        .filter(|site| is_hydrogen(site.element))
                        .map(|site| (residue.key.clone(), site.atom)),
                );
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
        let mut assumed = BTreeMap::new();
        for residue in &templated {
            raw::poll(self.cancellation)?;
            self.assign_template_charges(residue, &residues[&residue.key], &mut assumed)?;
        }
        self.report_assumed_protonation(assumed);
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

    /// Adds residues that took the usual pH 7 state to the structure's one
    /// report of them.
    fn report_assumed_protonation(&mut self, assumed: BTreeMap<String, usize>) {
        if assumed.is_empty() {
            return;
        }
        let existing = self.issues.iter_mut().find_map(|issue| match issue {
            raw::StructureInterpretIssue::ProtonationAssumed { components } => Some(components),
            _ => None,
        });
        let Some(components) = existing else {
            self.issues
                .push(raw::StructureInterpretIssue::ProtonationAssumed {
                    components: assumed.into_iter().collect(),
                });
            return;
        };
        let mut merged = components.drain(..).collect::<BTreeMap<_, _>>();
        for (component, residues) in assumed {
            *merged.entry(component).or_default() += residues;
        }
        components.extend(merged);
    }

    /// Charges the atoms of a residue that no template gives charges. Stated
    /// charges are kept, and the charges of the block's own component atom
    /// table apply. Otherwise charges follow from bonds. A source can state
    /// bond orders but no charges, as the component tables of wwPDB entry
    /// files do. A neutral atom with more bonds than its valence allows,
    /// such as the methylated N7 of 7-methylguanosine, then takes the charge
    /// of one that makes its valence valid, +1 before -1. A residue of one
    /// atom, or one holding a metal, that nothing gives charges is reported.
    fn assign_implied_charges(&mut self, residue: &ResidueAtoms) -> Result<(), FormatError> {
        let block = self
            .catalog
            .component_charges
            .get(&normalized(&residue.key.component_id));
        let stated = residue
            .sites
            .iter()
            .any(|site| site.stated_charge.is_some());
        if block.is_none()
            && !stated
            && (residue.sites.len() == 1
                || residue.sites.iter().any(|site| raw::is_metal(site.element)))
        {
            self.issues
                .push(raw::StructureInterpretIssue::TemplateMissing {
                    residue: residue.key.id(self.model_id),
                    component_id: residue.key.component_id.clone(),
                });
        }
        let mut charged_atoms = Vec::new();
        let mut unknown = Vec::new();
        for site in &residue.sites {
            if site.stated_charge.is_some() {
                continue;
            }
            if let Some(charge) = block.and_then(|atoms| atoms.get(&site.name)) {
                match charge {
                    Some(charge) => {
                        self.graph
                            .require_atom_mut(site.atom)
                            .map_err(interpret_error)?
                            .formal_charge = *charge;
                    }
                    None => unknown.push(site.name.clone()),
                }
                continue;
            }
            let atom = self
                .graph
                .require_atom(site.atom)
                .map_err(interpret_error)?
                .clone();
            if atom.formal_charge != 0
                || !matches!(
                    rdkit_valence_violation(self.graph, site.atom, &atom),
                    Some(ValenceViolation::ValenceExceeded { .. })
                )
            {
                continue;
            }
            let valid = |charge: i8| {
                let mut charged = atom.clone();
                charged.formal_charge = charge;
                rdkit_valence_violation(self.graph, site.atom, &charged).is_none()
            };
            let Some(charge) = [1, -1].into_iter().find(|&charge| valid(charge)) else {
                continue;
            };
            self.graph
                .require_atom_mut(site.atom)
                .map_err(interpret_error)?
                .formal_charge = charge;
            charged_atoms.push((site.name.clone(), charge));
            if charge == 1 {
                charged_atoms.extend(self.charge_terminal_partners(residue, site.atom, block)?);
            }
        }
        if !unknown.is_empty() {
            self.issues
                .push(raw::StructureInterpretIssue::TemplateChargeUnknown {
                    residue: residue.key.id(self.model_id),
                    component_id: residue.key.component_id.clone(),
                    atoms: unknown,
                });
        }
        if !charged_atoms.is_empty() {
            self.issues
                .push(raw::StructureInterpretIssue::ChargesFromBonds {
                    residue: residue.key.id(self.model_id),
                    component_id: residue.key.component_id.clone(),
                    charged_atoms,
                });
        }
        Ok(())
    }

    /// Charges -1 the terminal neighbours of an atom that took +1 from its
    /// bonds when, neutral, they would take a hydrogen: the oxygens of a
    /// nitro group, a nitrate, or an N-oxide, and the end of an azide. The
    /// bonds drawn as N=O and N-O then hold the charge-separated form, not
    /// an N-OH.
    fn charge_terminal_partners(
        &mut self,
        residue: &ResidueAtoms,
        cation: AtomId,
        block: Option<&BTreeMap<String, Option<i8>>>,
    ) -> Result<Vec<(String, i8)>, FormatError> {
        let mut charged = Vec::new();
        let neighbors = self
            .graph
            .require_neighbors(cation)
            .map_err(interpret_error)?
            .collect::<Vec<_>>();
        for neighbor in neighbors {
            let Some(site) = residue.site(neighbor) else {
                continue;
            };
            let blocked = block.is_some_and(|atoms| atoms.contains_key(&site.name));
            if site.stated_charge.is_some()
                || blocked
                || !matches!(site.element.atomic_number(), 7 | 8 | 16)
                || self
                    .graph
                    .require_neighbors(neighbor)
                    .map_err(interpret_error)?
                    .count()
                    != 1
            {
                continue;
            }
            let atom = self
                .graph
                .require_atom(neighbor)
                .map_err(interpret_error)?
                .clone();
            if atom.formal_charge != 0
                || rdkit_inferred_hydrogen_count(self.graph, neighbor, &atom) == 0
            {
                continue;
            }
            let mut anion = atom;
            anion.formal_charge = -1;
            if rdkit_valence_violation(self.graph, neighbor, &anion).is_some()
                || rdkit_inferred_hydrogen_count(self.graph, neighbor, &anion) != 0
            {
                continue;
            }
            self.graph
                .require_atom_mut(neighbor)
                .map_err(interpret_error)?
                .formal_charge = -1;
            charged.push((site.name.clone(), -1));
        }
        Ok(charged)
    }

    /// Whether a declared connection already bonds two atoms of the residue.
    fn has_stated_internal_bond(&self, residue: &ResidueAtoms) -> Result<bool, FormatError> {
        for site in &residue.sites {
            if self
                .graph
                .require_neighbors(site.atom)
                .map_err(interpret_error)?
                .any(|neighbor| residue.site(neighbor).is_some())
            {
                return Ok(true);
            }
        }
        Ok(false)
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

    /// Adds listed bonds between atoms not yet bonded; an existing bond keeps
    /// its order, since connectivity-only records state none reliably.
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
        template: Option<&(&Component, Vec<&Component>)>,
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
        let templated = template.and_then(|(component, variants)| {
            std::iter::once(*component)
                .chain(variants.iter().copied())
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
    /// templates when the source does not describe the entity.
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
            // The flag states that at least one monomer-to-monomer linkage
            // differs from the linkage implied by the polymer type, but does
            // not identify that pair. Remain conservative rather than
            // inventing a standard bond.
            if entity
                .as_ref()
                .is_some_and(|entity| self.catalog.nonstandard_polymer_linkage.contains(entity))
            {
                continue;
            }
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
                let Some((&left_atom, &right_atom)) = link.and_then(|link| {
                    link.atoms(|name| left.atoms.get(name), |name| right.atoms.get(name))
                }) else {
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
                        SourceBondOrder::Localized(BondOrder::Single),
                        raw::BondEvidence::PolymerLink,
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

    /// Bonds each unbonded hydrogen that its template or listed bonds do not
    /// cover. A hydrogen that the template's variants name, or a terminal
    /// hydrogen by convention (see [`terminal_hydrogen_anchors`]), bonds to
    /// the atom its name places it on when that atom lies within bonding
    /// distance; any other bonds to the nearest heavy atom of its residue
    /// within bonding distance, which is reported.
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
            if let Some((anchor, evidence)) = self.named_anchor(residue, site) {
                self.add_bond(
                    hydrogen,
                    anchor,
                    SourceBondOrder::Localized(BondOrder::Single),
                    evidence,
                )?;
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

    /// The atom a hydrogen's name places it on, with the evidence: the
    /// atom that every variant naming the hydrogen bonds it to, or a
    /// terminal hydrogen's atom. `None` when the name places it nowhere or
    /// that atom is absent or beyond bonding distance.
    fn named_anchor(
        &self,
        residue: &ResidueAtoms,
        hydrogen: &ResidueSite,
    ) -> Option<(AtomId, raw::BondEvidence)> {
        let (component, variants) = self
            .catalog
            .template(&normalized(&residue.key.component_id))?;
        let mut named = variants.iter().filter_map(|variant| {
            let (index, atom) = variant.atom(&hydrogen.name)?;
            if !is_hydrogen(atom.element()) {
                return None;
            }
            let bond = variant
                .bonds()
                .iter()
                .find(|bond| bond.atoms().contains(&index))?;
            let [a, b] = bond.atoms();
            let anchor = variant.atoms()[if a == index { b } else { a }].name();
            Some((anchor, variant.variant()))
        });
        let (anchors, evidence) = match named.next() {
            Some((anchor, variant)) => {
                if named.any(|(other, _)| other != anchor) {
                    return None;
                }
                (
                    vec![anchor],
                    raw::BondEvidence::Template {
                        component_id: component.id().to_owned(),
                        variant: variant.map(str::to_owned),
                    },
                )
            }
            None => (
                terminal_hydrogen_anchors(component.linking(), &hydrogen.name).to_vec(),
                raw::BondEvidence::HydrogenName,
            ),
        };
        let anchor = anchors
            .iter()
            .find_map(|name| residue.atoms.get(*name))
            .and_then(|&atom| residue.site(atom))?;
        within_bonding_distance(
            (hydrogen.element, hydrogen.point),
            (anchor.element, anchor.point),
        )
        .then_some((anchor.atom, evidence))
    }

    /// Sets the formal charges of a templated residue from the protonation
    /// state its hydrogens, links, and metal coordination fit (see
    /// [`protonation`]), checks stated charges against that state, and
    /// reports what the template alone does not give: missing atoms and
    /// hydrogens, unknown template charges, and charges the hydrogens or
    /// coordination change. `assumed` counts the residues, by component,
    /// that take the usual pH 7 state for want of hydrogens.
    fn assign_template_charges(
        &mut self,
        templated: &TemplatedResidue<'_>,
        residue: &ResidueAtoms,
        assumed: &mut BTreeMap<String, usize>,
    ) -> Result<(), FormatError> {
        let component = templated.component;
        let named = templated
            .matched
            .iter()
            .map(|(&index, &atom)| (component.atoms()[index].name(), atom))
            .collect::<Vec<_>>();
        let mut heavy = BTreeMap::new();
        for &(name, atom) in &named {
            let site = residue.site(atom).expect("matched atom is a residue site");
            if !is_hydrogen(site.element) {
                heavy.insert(name, self.atom_evidence(atom, residue)?);
            }
        }
        let evidence = protonation::ResidueEvidence {
            heavy,
            has_hydrogens: residue.sites.iter().any(|site| is_hydrogen(site.element)),
            usual_state: templated.usual_state,
        };
        let candidates = std::iter::once(component).chain(templated.variants.iter().copied());
        let chosen = protonation::choose(candidates, &evidence).ok_or_else(|| {
            raw::structure_error(
                DiagnosticKind::Inconsistent,
                residue.sites.first().map(|site| site.line),
                format!(
                    "The hydrogens of a {:?} residue match neither its template nor a protonation variant.",
                    component.id()
                ),
            )
        })?;
        let candidate = chosen.candidate;
        if candidate.variant().is_some() {
            self.apply_variant_bonds(candidate, &named)?;
        }
        let atom_of = |name: &str| {
            named
                .iter()
                .find(|(candidate, _)| *candidate == name)
                .map(|&(_, atom)| atom)
                .expect("states cover matched atoms")
        };
        let residue_id = templated.key.id(self.model_id);
        let component_id = component.id().to_owned();
        let (mut from_hydrogens, mut missing, mut coordinated, mut unknown) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for state in &chosen.atoms {
            let atom = atom_of(&state.name);
            let site = residue.site(atom).expect("matched atom is a residue site");
            let mut value = state.charge;
            let mut inferred = site.stated_charge.is_none();
            if let Some(stated) = site.stated_charge.filter(|&stated| stated != state.charge) {
                if self.holds(atom, stated, state.hydrogens)? {
                    value = stated;
                    self.issues
                        .push(raw::StructureInterpretIssue::StatedChargeKept {
                            residue: residue_id.clone(),
                            component_id: component_id.clone(),
                            atom_name: site.name.clone(),
                            stated,
                            derived: state.charge,
                        });
                } else {
                    inferred = true;
                    self.issues
                        .push(raw::StructureInterpretIssue::StatedChargeRejected {
                            residue: residue_id.clone(),
                            component_id: component_id.clone(),
                            atom_name: site.name.clone(),
                            stated,
                            kept: state.charge,
                        });
                }
            }
            self.graph
                .require_atom_mut(atom)
                .map_err(interpret_error)?
                .formal_charge = value;
            if inferred && state.from_hydrogens != 0 {
                from_hydrogens.push((site.name.clone(), value));
            }
            if state.missing > 0 {
                missing.push(site.name.clone());
            }
            if state.coordinated > 0 {
                coordinated.push(site.name.clone());
            }
            if state.charge_unknown && site.stated_charge.is_none() {
                unknown.push(site.name.clone());
            }
        }
        if !chosen.without_evidence && (candidate.variant().is_some() || !from_hydrogens.is_empty())
        {
            self.issues
                .push(raw::StructureInterpretIssue::ProtonationFromHydrogens {
                    residue: residue_id.clone(),
                    component_id: component_id.clone(),
                    variant: candidate.variant().map(str::to_owned),
                    charged_atoms: from_hydrogens,
                });
        }
        // A residue with no heteroatom, such as a CA trace, assumes nothing.
        let heteroatoms = residue
            .sites
            .iter()
            .any(|site| !is_hydrogen(site.element) && site.element.atomic_number() != 6);
        if chosen.without_evidence && templated.usual_state && heteroatoms {
            *assumed.entry(component_id.clone()).or_default() += 1;
        }
        if !missing.is_empty() {
            self.issues
                .push(raw::StructureInterpretIssue::HydrogensMissing {
                    residue: residue_id.clone(),
                    component_id: component_id.clone(),
                    atoms: missing,
                });
        }
        if !coordinated.is_empty() {
            self.issues.push(
                raw::StructureInterpretIssue::HydrogenReplacedByCoordination {
                    residue: residue_id.clone(),
                    component_id: component_id.clone(),
                    atoms: coordinated,
                },
            );
        }
        if !unknown.is_empty() {
            self.issues
                .push(raw::StructureInterpretIssue::TemplateChargeUnknown {
                    residue: residue_id,
                    component_id,
                    atoms: unknown,
                });
        }
        self.report_missing_atoms(templated, &evidence, &named)
    }

    /// Reports the template heavy atoms a residue lacks. A leaving atom
    /// that nothing replaces, such as the `OXT` of a residue before a gap,
    /// leaves its bonded atom without the hydrogen that valence would add,
    /// so the residue does not become another molecule. A nucleotide
    /// without its whole 5' phosphate is a 5'-hydroxyl terminus, not a
    /// residue with missing atoms.
    fn report_missing_atoms(
        &mut self,
        templated: &TemplatedResidue<'_>,
        evidence: &protonation::ResidueEvidence<'_>,
        named: &[(&str, AtomId)],
    ) -> Result<(), FormatError> {
        let component = templated.component;
        let has = |name: &str| evidence.heavy.contains_key(name);
        let five_prime_hydroxyl = component.linking() == ComponentLinking::Nucleotide
            && component
                .atom("P")
                .is_some_and(|(_, atom)| !has(atom.name()));
        let mut missing = Vec::new();
        for (index, atom) in component.atoms().iter().enumerate() {
            if is_hydrogen(atom.element()) || templated.matched.contains_key(&index) {
                continue;
            }
            if five_prime_hydroxyl && FIVE_PRIME_PHOSPHATE.contains(&atom.name()) {
                continue;
            }
            if !atom.leaving() {
                missing.push(atom.name().to_owned());
                continue;
            }
            for bond in component.bonds() {
                let [a, b] = bond.atoms();
                let anchor = match (a == index, b == index) {
                    (true, _) => b,
                    (_, true) => a,
                    _ => continue,
                };
                let anchor_name = component.atoms()[anchor].name();
                let Some(anchor_evidence) = evidence.heavy.get(anchor_name) else {
                    continue;
                };
                if anchor_evidence.links > 0 {
                    continue;
                }
                let anchor_atom = named
                    .iter()
                    .find(|(name, _)| *name == anchor_name)
                    .map(|&(_, atom)| atom)
                    .expect("evidence covers matched atoms");
                let fixed = fixed_hydrogens(
                    component,
                    anchor_name,
                    |atom| has(atom.name()),
                    false,
                    anchor_evidence.hydrogens,
                )
                .expect("the anchor lacks a leaving atom");
                self.graph
                    .require_atom_mut(anchor_atom)
                    .map_err(interpret_error)?
                    .hydrogens = ImplicitHydrogens::Fixed(fixed);
                self.issues
                    .push(raw::StructureInterpretIssue::TerminalAtomMissing {
                        residue: templated.key.id(self.model_id),
                        component_id: component.id().to_owned(),
                        atom_name: atom.name().to_owned(),
                        anchor: anchor_name.to_owned(),
                    });
            }
        }
        if !missing.is_empty() {
            self.issues
                .push(raw::StructureInterpretIssue::AtomsMissing {
                    residue: templated.key.id(self.model_id),
                    component_id: component.id().to_owned(),
                    atoms: missing,
                });
        }
        Ok(())
    }

    /// What the model shows about one matched heavy atom.
    fn atom_evidence(
        &self,
        atom: AtomId,
        residue: &ResidueAtoms,
    ) -> Result<protonation::AtomEvidence, FormatError> {
        let hydrogens = self
            .graph
            .require_neighbors(atom)
            .map_err(interpret_error)?
            .filter(|&neighbor| {
                self.graph
                    .atom(neighbor)
                    .is_some_and(|neighbor| is_hydrogen(neighbor.element))
            })
            .count();
        Ok(protonation::AtomEvidence {
            hydrogens,
            links: self.links_outside(atom, residue)?,
            coordination: self.coordination.get(&atom).copied().unwrap_or(0),
        })
    }

    /// Whether an atom's bonds and hydrogens can hold a stated charge: its
    /// valence is valid, and it needs no more hydrogens than it has or the
    /// protonation state gives it (`state_hydrogens`). A +1 on a
    /// carboxylate oxygen would need two hydrogens and is not held.
    fn holds(&self, atom: AtomId, stated: i8, state_hydrogens: usize) -> Result<bool, FormatError> {
        let mut charged = self
            .graph
            .require_atom(atom)
            .map_err(interpret_error)?
            .clone();
        charged.formal_charge = stated;
        if rdkit_valence_violation(self.graph, atom, &charged).is_some() {
            return Ok(false);
        }
        let explicit = self
            .graph
            .require_neighbors(atom)
            .map_err(interpret_error)?
            .filter(|&neighbor| {
                self.graph
                    .atom(neighbor)
                    .is_some_and(|neighbor| is_hydrogen(neighbor.element))
            })
            .count();
        let inferred = usize::from(rdkit_inferred_hydrogen_count(self.graph, atom, &charged));
        Ok(explicit + inferred <= explicit.max(state_hydrogens))
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

    /// Replaces template bond orders with a variant's.
    fn apply_variant_bonds(
        &mut self,
        variant: &Component,
        named: &[(&str, AtomId)],
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
