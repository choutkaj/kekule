//! Resolves the bonds a macromolecular source does not state, from ranked
//! evidence: the source's own component bond table or connection records,
//! component templates, polymer and branch linkage, and, for residues nothing
//! else covers, covalent radii.
use crate::diagnostics::{DiagnosticKind, FormatError};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use crate::chemistry::localize_source_aromatic_bonds;
use crate::components::{
    Component, ComponentBond as TemplateBond, ComponentBondOrder as TemplateOrder,
    ComponentLibrary, ComponentLinking,
};
use crate::geometry::Point3;
use crate::molecule::{AtomId, BondOrder, Element, Molecule, MoleculeEditor};

use super::mmcif_category::MmcifCategory;
use super::structure_interpret as raw;
use super::MmcifBlock;

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

/// Resolves the bonds of one staged molecule group in place. Bonds already
/// present, from connection records, must agree with every other source.
pub(crate) fn complete_editor_connectivity<'a>(
    catalog: &ConnectivityCatalog,
    editor: &mut MoleculeEditor,
    atoms: impl IntoIterator<Item = StagedAtomProvenance<'a>>,
    listed: &[ListedBond],
    model_id: &str,
    evidence: &mut BondEvidenceMap,
    issues: &mut Vec<raw::StructureInterpretIssue>,
) -> Result<(), FormatError> {
    let mut source_aromatic_bonds = BTreeSet::new();
    let mut resolver = Resolver {
        catalog,
        graph: editor.working_mut(),
        source_aromatic_bonds: &mut source_aromatic_bonds,
        evidence,
        issues,
        model_id,
    };
    resolver.resolve(residue_atoms(atoms), listed)?;
    localize_source_aromatic_bonds(editor.working_mut(), &source_aromatic_bonds)
        .map_err(interpret_error)
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

    /// Metal coordination is covalent but does not count toward the valence
    /// of the ligand atom, so it is a dative bond whatever its stated order.
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
}

fn residue_atoms<'a>(
    atoms: impl IntoIterator<Item = StagedAtomProvenance<'a>>,
) -> BTreeMap<ResidueKey, ResidueAtoms> {
    let mut residues = BTreeMap::new();
    for atom in atoms {
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
        let residue = residues.entry(key.clone()).or_insert_with(|| ResidueAtoms {
            key,
            atoms: BTreeMap::new(),
            sites: Vec::new(),
        });
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
    residues
}

const COVALENT_BOND_TOLERANCE_ANGSTROM: f64 = 0.45;
const MIN_COVALENT_BOND_DISTANCE_SQUARED: f64 = 0.16;
const FALLBACK_COVALENT_RADIUS_ANGSTROM: f64 = 0.77;

/// Whether two atoms lie within covalent bonding distance: the sum of their
/// covalent radii plus a tolerance, and not implausibly close.
pub(crate) fn within_bonding_distance(
    (left, left_point): (Element, Point3),
    (right, right_point): (Element, Point3),
) -> bool {
    let radius = |element: Element| {
        element
            .covalent_radius_angstrom()
            .unwrap_or(FALLBACK_COVALENT_RADIUS_ANGSTROM)
    };
    let cutoff = radius(left) + radius(right) + COVALENT_BOND_TOLERANCE_ANGSTROM;
    let distance_squared = (left_point.x - right_point.x).powi(2)
        + (left_point.y - right_point.y).powi(2)
        + (left_point.z - right_point.z).powi(2);
    distance_squared > MIN_COVALENT_BOND_DISTANCE_SQUARED && distance_squared <= cutoff * cutoff
}

/// The coordination flag `component` gives the bond between the two named
/// atoms of a component-table row, if it lists that bond.
fn template_coordination(component: &Component, bond: &ComponentBond) -> Option<bool> {
    let (left, _) = component.atom(&bond.atom_1)?;
    let (right, _) = component.atom(&bond.atom_2)?;
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
    graph: &'g mut Molecule,
    source_aromatic_bonds: &'g mut BTreeSet<crate::molecule::BondId>,
    evidence: &'g mut BondEvidenceMap,
    issues: &'g mut Vec<raw::StructureInterpretIssue>,
    model_id: &'r str,
}

impl<'r> Resolver<'r, '_> {
    fn resolve(
        &mut self,
        residues: BTreeMap<ResidueKey, ResidueAtoms>,
        listed: &[ListedBond],
    ) -> Result<(), FormatError> {
        let mut residue_of = BTreeMap::new();
        for residue in residues.values() {
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
            let listed = listed_in.remove(&residue.key).unwrap_or_default();
            let component_id = normalized(&residue.key.component_id);
            let template = self.catalog.template(&component_id);
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
        }
        self.add_listed_bonds(&listed_across)?;
        self.apply_branch_links(&residues)?;
        self.apply_polymer_links(&residues)?;
        self.attach_stray_hydrogens(&residues, stray_hydrogens)?;
        for residue in &templated {
            self.assign_template_charges(residue, &residues[&residue.key])?;
        }
        self.report_fragmented_residues(&residues);
        Ok(())
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
        let templated = template.and_then(|(component, variants)| {
            std::iter::once(*component)
                .chain(variants.iter().copied())
                .find_map(|component| template_coordination(component, bond))
        });
        templated.unwrap_or_else(|| {
            let [left, right] = atoms.map(|atom| {
                residue
                    .site(atom)
                    .expect("component-table atoms are residue sites")
                    .element
            });
            raw::is_metal(left) != raw::is_metal(right)
        })
    }

    /// Bonds every pair of a residue's atoms within covalent bonding
    /// distance, when the caller allows it.
    fn infer_residue_bonds(&mut self, residue: &ResidueAtoms) -> Result<(), FormatError> {
        if !self.catalog.infer_bonds_from_distance {
            return Err(raw::structure_error(
                DiagnosticKind::Unsupported,
                residue.sites.first().map(|site| site.line),
                format!(
                    "component `{}` has no template or bond table; supply a template or enable distance inference",
                    residue.key.component_id
                ),
            ));
        }
        let mut bonds = 0;
        for (index, left) in residue.sites.iter().enumerate() {
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
            for pair in chain.windows(2) {
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
                let peptide = || left.atoms.get("C").zip(right.atoms.get("N"));
                let nucleotide = || {
                    left.atoms
                        .get("O3'")
                        .or_else(|| left.atoms.get("O3*"))
                        .zip(right.atoms.get("P"))
                };
                let pair = match declared {
                    Some(kind) if is_polypeptide(kind) => peptide(),
                    Some(kind) if is_nucleic_acid(kind) => nucleotide(),
                    Some(_) => None,
                    // Without a declared polymer type, the residues' templates
                    // name the linkage. A residue without a polymer template
                    // defers to its neighbour, and two such residues link by
                    // atom names; bonding distance decides either way.
                    None => {
                        let polymer = |linking| {
                            matches!(
                                linking,
                                ComponentLinking::Peptide | ComponentLinking::Nucleotide
                            )
                            .then_some(linking)
                        };
                        match (
                            polymer(self.catalog.linking(&left.key.component_id)),
                            polymer(self.catalog.linking(&right.key.component_id)),
                        ) {
                            (Some(first), Some(second)) if first != second => None,
                            (Some(ComponentLinking::Peptide), _)
                            | (_, Some(ComponentLinking::Peptide)) => peptide(),
                            (Some(_), _) | (_, Some(_)) => nucleotide(),
                            (None, None) => peptide().or_else(nucleotide),
                        }
                    }
                };
                let Some((&left_atom, &right_atom)) = pair else {
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
    /// cover to the nearest heavy atom of its residue within bonding distance.
    fn attach_stray_hydrogens(
        &mut self,
        residues: &BTreeMap<ResidueKey, ResidueAtoms>,
        stray_hydrogens: Vec<(ResidueKey, AtomId)>,
    ) -> Result<(), FormatError> {
        for (key, hydrogen) in stray_hydrogens {
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

    /// Sets formal charges of a templated residue. Stated charges win. Without
    /// model hydrogens, the template's charges apply. With them, the template
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
        let candidates =
            std::iter::once(templated.component).chain(templated.variants.iter().copied());
        let mut evaluated = Vec::new();
        for candidate in candidates {
            if let Some(charges) =
                self.candidate_charges(candidate, &named, residue, has_hydrogens)?
            {
                evaluated.push((candidate, charges));
            }
        }
        let chosen = evaluated
            .iter()
            .find(|(_, charges)| charges.iter().all(|charge| charge.delta == 0))
            .or_else(|| evaluated.first())
            .ok_or_else(|| {
                raw::structure_error(
                    DiagnosticKind::Inconsistent,
                    residue.sites.first().map(|site| site.line),
                    format!(
                        "hydrogens of a `{}` residue match neither its template nor a protonation variant",
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
        if candidate.variant().is_some() {
            self.apply_variant_bonds(candidate, &named)?;
        }
        if candidate.variant().is_some() || !charged_atoms.is_empty() {
            self.issues
                .push(raw::StructureInterpretIssue::ProtonationFromHydrogens {
                    residue: templated.key.id(self.model_id),
                    component_id: templated.component.id().to_owned(),
                    variant: candidate.variant().map(str::to_owned),
                    charged_atoms,
                });
        }
        Ok(())
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
        let mut charges = Vec::new();
        for &(name, atom) in named {
            let Some((index, template_atom)) = candidate.atom(name) else {
                // Variants only remove hydrogens; a model hydrogen the
                // variant lacks shows up on its heavy atom below.
                if is_hydrogen(residue.site(atom).expect("matched site").element) {
                    continue;
                }
                return Ok(None);
            };
            let base = template_atom.formal_charge().unwrap_or(0);
            let element = template_atom.element();
            let mut delta = 0i8;
            if has_hydrogens && !is_hydrogen(element) && !is_carbon(element) {
                let expected =
                    expected_hydrogens(candidate, index, self.links_outside(atom, residue)?);
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
                    let difference = present as isize - expected as isize;
                    // A ring nitrogen takes or loses one proton whatever its
                    // bond orders, as a protonated cytosine N3 or a
                    // deprotonated imidazole nitrogen does.
                    let ring_nitrogen = element.atomic_number() == 7
                        && difference.abs() == 1
                        && in_ring(candidate, index);
                    if !only_single && !ring_nitrogen {
                        return Ok(None);
                    }
                    delta = i8::try_from(difference).map_err(|_| {
                        chemistry_error("hydrogen count difference exceeds charge range")
                    })?;
                }
            }
            charges.push(TemplateCharge {
                atom,
                value: base.saturating_add(delta),
                delta,
            });
        }
        Ok(Some(charges))
    }

    /// Whether a residue atom is bonded to an atom of another residue.
    fn links_outside(&self, atom: AtomId, residue: &ResidueAtoms) -> Result<bool, FormatError> {
        Ok(self
            .graph
            .require_neighbors(atom)
            .map_err(interpret_error)?
            .any(|neighbor| residue.site(neighbor).is_none()))
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
        let mut component_of = BTreeMap::new();
        for (index, component) in components.iter().enumerate() {
            for atom in component {
                component_of.insert(*atom, index);
            }
        }
        for residue in residues.values() {
            let fragments = residue
                .sites
                .iter()
                .filter_map(|site| component_of.get(&site.atom))
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
                        "atom `{}` of component `{}` is {} in the model but {} in its template",
                        site.name,
                        component.id(),
                        site.element.symbol(),
                        template_atom.element().symbol()
                    ));
                }
                if matched.insert(index, site.atom).is_some() {
                    return failure(format!(
                        "two atoms of one `{}` residue match template atom `{}`",
                        component.id(),
                        template_atom.name()
                    ));
                }
            }
            None if is_hydrogen(site.element) => strays.push(site.atom),
            None => {
                return failure(format!(
                    "atom `{}` is not in component `{}`",
                    site.name,
                    component.id()
                ))
            }
        }
    }
    Ok(Some((matched, strays)))
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

/// The hydrogens a template atom carries, less leaving hydrogens when the
/// atom links to another residue.
fn expected_hydrogens(component: &Component, index: usize, linked: bool) -> usize {
    component
        .bonds()
        .iter()
        .filter_map(|bond| {
            let [a, b] = bond.atoms();
            match (a == index, b == index) {
                (true, _) => Some(b),
                (_, true) => Some(a),
                _ => None,
            }
        })
        .filter(|&other| {
            let atom = &component.atoms()[other];
            is_hydrogen(atom.element()) && !(linked && atom.leaving())
        })
        .count()
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
    graph: &mut Molecule,
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
            SourceBondOrder::Localized(order) => format!("{order:?}"),
            SourceBondOrder::Aromatic => "Aromatic".to_owned(),
        };
        return Err(chemistry_error(format!(
            "conflicting authoritative bond evidence for one atom pair: {existing_order:?} versus {source_order}"
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
            format!("unsupported component bond order `{other}`"),
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
            format!("missing required {tag}"),
        )
    })?;
    value.optional_text().ok_or_else(|| {
        row_error(
            table,
            row,
            DiagnosticKind::InvalidValue,
            format!("missing required {tag}"),
        )
    })
}

fn required_i32(table: &MmcifCategory<'_>, row: usize, tag: &str) -> Result<i32, FormatError> {
    required(table, row, tag)?.parse::<i32>().map_err(|_| {
        row_error(
            table,
            row,
            DiagnosticKind::InvalidValue,
            format!("invalid integer {tag}"),
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
    chemistry_error("molecule construction failed").with_source(source)
}

fn chemistry_error(message: impl Into<String>) -> FormatError {
    raw::structure_error(DiagnosticKind::Chemistry, None, message)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::MoleculeInstanceId;

    const PEPTIDE: &str = r#"
data_peptide
loop_
_entity.id
_entity.type
1 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.nstd_linkage
1 'polypeptide(L)' no
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.mon_id
A 1 1 GLY
A 1 2 GLY
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
GLY N CA sing
GLY CA C sing
GLY C O doub
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
ATOM 1 N N GLY A 1 1 0.0 0.0 0.0
ATOM 2 C CA GLY A 1 1 1.4 0.0 0.0
ATOM 3 C C GLY A 1 1 2.8 0.0 0.0
ATOM 4 O O GLY A 1 1 3.8 0.0 0.0
ATOM 5 N N GLY A 1 2 4.1 0.0 0.0
ATOM 6 C CA GLY A 1 2 5.5 0.0 0.0
ATOM 7 C C GLY A 1 2 6.9 0.0 0.0
ATOM 8 O O GLY A 1 2 7.9 0.0 0.0
"#;

    const BRANCH: &str = r#"
data_branch
loop_
_entity.id
_entity.type
2 branched
loop_
_struct_asym.id
_struct_asym.entity_id
B 2
loop_
_pdbx_entity_branch_list.entity_id
_pdbx_entity_branch_list.comp_id
_pdbx_entity_branch_list.num
_pdbx_entity_branch_list.hetero
2 NAG 1 n
2 GAL 2 n
loop_
_pdbx_branch_scheme.asym_id
_pdbx_branch_scheme.entity_id
_pdbx_branch_scheme.mon_id
_pdbx_branch_scheme.num
_pdbx_branch_scheme.pdb_seq_num
B 2 NAG 1 10
B 2 GAL 2 11
loop_
_pdbx_entity_branch_link.link_id
_pdbx_entity_branch_link.entity_id
_pdbx_entity_branch_link.entity_branch_list_num_1
_pdbx_entity_branch_link.comp_id_1
_pdbx_entity_branch_link.atom_id_1
_pdbx_entity_branch_link.entity_branch_list_num_2
_pdbx_entity_branch_link.comp_id_2
_pdbx_entity_branch_link.atom_id_2
_pdbx_entity_branch_link.value_order
1 2 1 NAG O4 2 GAL C1 sing
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
NAG C4 O4 sing
GAL C1 O5 sing
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C4 NAG B 2 . 10 0.0 0.0 0.0
HETATM 2 O O4 NAG B 2 . 10 1.4 0.0 0.0
HETATM 3 C C1 GAL B 2 . 11 2.8 0.0 0.0
HETATM 4 O O5 GAL B 2 . 11 4.2 0.0 0.0
"#;

    #[test]
    fn singleton_connectivity_categories_are_syntax_independent() {
        for (input, bonds) in [(PEPTIDE, 7), (BRANCH, 3)] {
            let scalar = crate::tests::mmcif_syntax::singleton_loops_as_scalars(input);
            let a = crate::mmcif::parse_str(input, Default::default())
                .unwrap()
                .interpret(Default::default())
                .unwrap();
            let b = crate::mmcif::parse_str(&scalar, Default::default())
                .unwrap()
                .interpret(Default::default())
                .unwrap();
            assert_eq!(b.topology().bond_count(), bonds);
            assert_eq!(b.topology().instance_count(), 1);
            assert!(a.topology().same_layout(b.topology()));
            assert_eq!(a.model().positions(), b.model().positions());
        }
    }

    #[test]
    fn singleton_component_bond_and_branch_site_are_read_as_scalars() {
        let source = "data_x\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nGLY C O doub\nloop_\n_pdbx_branch_scheme.entity_id\n_pdbx_branch_scheme.asym_id\n_pdbx_branch_scheme.num\n_pdbx_branch_scheme.mon_id\n_pdbx_branch_scheme.auth_seq_num\n2 B 1 NAG 10\n";
        for input in [
            source.to_owned(),
            crate::tests::mmcif_syntax::singleton_loops_as_scalars(source),
        ] {
            let doc = crate::mmcif::parse_str(&input, Default::default()).unwrap();
            let catalog = ConnectivityCatalog::from_block(
                &doc.blocks()[0],
                &raw::ConnectivityOptions::default(),
            )
            .unwrap();
            let bond = &catalog.component_bonds["GLY"][0];
            assert_eq!((&*bond.atom_1, &*bond.atom_2), ("C", "O"));
            assert!(matches!(
                bond.order,
                SourceBondOrder::Localized(BondOrder::Double)
            ));
            assert_eq!(catalog.branch_sites.len(), 1);
            let site = &catalog.branch_sites[0];
            assert_eq!(
                (
                    &*site.entity_id,
                    &*site.asym_id,
                    site.number,
                    &*site.component_id
                ),
                ("2", "B", 1, "NAG")
            );
            assert_eq!(site.author_sequence_id.as_deref(), Some("10"));
        }
    }

    #[test]
    fn component_templates_and_polymer_links_complete_peptide_graph() {
        let document = super::super::mmcif_document::parse_mmcif_str(
            PEPTIDE,
            super::super::mmcif_document::MmcifParseOptions::default(),
        )
        .expect("parse peptide");
        let interpretation =
            raw::interpret_mmcif(&document, raw::StructureInterpretOptions::default())
                .expect("interpret peptide");
        let graph = interpretation
            .topology()
            .instance(MoleculeInstanceId::new(0))
            .expect("peptide instance")
            .molecule();
        assert_eq!(graph.bond_count(), 7);
        assert_eq!(
            graph
                .bonds()
                .filter(|(_, bond)| bond.order == BondOrder::Single)
                .count(),
            5
        );
        assert_eq!(
            graph
                .bonds()
                .filter(|(_, bond)| bond.order == BondOrder::Double)
                .count(),
            2
        );
        assert!(graph.is_connected());
    }

    #[test]
    fn sequence_gap_is_not_bridged_by_a_fake_polymer_bond() {
        let input = PEPTIDE
            .replace("A 1 2 GLY", "A 1 3 GLY")
            .replace("GLY A 1 2 4.1", "GLY A 1 3 4.1")
            .replace("GLY A 1 2 5.5", "GLY A 1 3 5.5")
            .replace("GLY A 1 2 6.9", "GLY A 1 3 6.9")
            .replace("GLY A 1 2 7.9", "GLY A 1 3 7.9");
        let document = super::super::mmcif_document::parse_mmcif_str(
            &input,
            super::super::mmcif_document::MmcifParseOptions::default(),
        )
        .expect("parse gapped peptide");
        let interpretation =
            raw::interpret_mmcif(&document, raw::StructureInterpretOptions::default())
                .expect("interpret gapped peptide");
        assert_eq!(interpretation.topology().instance_count(), 2);
        assert!(interpretation
            .topology()
            .instance_records()
            .all(|(instance, _)| interpretation
                .topology()
                .instance(instance)
                .expect("partitioned peptide instance")
                .molecule()
                .is_connected()));
    }

    #[test]
    fn branch_link_connects_branched_macromolecule() {
        let document = super::super::mmcif_document::parse_mmcif_str(
            BRANCH,
            super::super::mmcif_document::MmcifParseOptions::default(),
        )
        .expect("parse branch");
        let interpretation =
            raw::interpret_mmcif(&document, raw::StructureInterpretOptions::default())
                .expect("interpret branch");
        let graph = interpretation
            .topology()
            .instance(MoleculeInstanceId::new(0))
            .expect("branch instance")
            .molecule();
        assert_eq!(graph.bond_count(), 3);
        assert!(graph
            .bonds()
            .all(|(_, bond)| bond.order == BondOrder::Single));
        assert!(graph.is_connected());
    }

    /// Two copies, B and C, of [`BRANCH`]'s NAG-GAL entity, which a covalent
    /// `_struct_conn` from B's `O5` to C's `C4` joins into one molecule
    /// group. `rows` are the `_atom_site` rows of [`BRANCH`]'s layout.
    fn two_branch_copies(rows: &str) -> String {
        let header = BRANCH
            .split("loop_\n_atom_site.group_PDB")
            .next()
            .unwrap()
            .replace("B 2\n", "B 2\nC 2\n")
            .replace(
                "B 2 GAL 2 11\n",
                "B 2 GAL 2 11\nC 2 NAG 1 20\nC 2 GAL 2 21\n",
            );
        let atom_site = BRANCH.split("loop_\n_atom_site.group_PDB").nth(1).unwrap();
        let columns = atom_site.split("HETATM").next().unwrap();
        format!(
            "{header}loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n\
             _struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_comp_id\n\
             _struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_label_seq_id\n\
             _struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_comp_id\n\
             _struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_label_seq_id\n\
             covale1 covale B GAL O5 . C NAG C4 .\n\
             loop_\n_atom_site.group_PDB{columns}{rows}"
        )
    }

    #[test]
    fn branch_links_bond_each_copy_within_bonding_distance() {
        let copy_b = "HETATM 1 C C4 NAG B 2 . 10 0.0 0.0 0.0\n\
                      HETATM 2 O O4 NAG B 2 . 10 1.4 0.0 0.0\n";
        let copy_b_gal = "HETATM 3 C C1 GAL B 2 . 11 2.8 0.0 0.0\n\
                          HETATM 4 O O5 GAL B 2 . 11 4.2 0.0 0.0\n";
        let copy_c = "HETATM 5 C C4 NAG C 2 . 20 5.6 0.0 0.0\n\
                      HETATM 6 O O4 NAG C 2 . 20 7.0 0.0 0.0\n";
        let copy_c_gal = "HETATM 7 C C1 GAL C 2 . 21 8.4 0.0 0.0\n\
                          HETATM 8 O O5 GAL C 2 . 21 9.8 0.0 0.0\n";
        // C's galactose 5 Å off its link partner.
        let copy_c_gal_far = "HETATM 7 C C1 GAL C 2 . 21 8.4 5.0 0.0\n\
                              HETATM 8 O O5 GAL C 2 . 21 9.8 5.0 0.0\n";
        let residue = |asym: &str, sequence: &str| SourceResidueId {
            model_id: "1".to_owned(),
            asym_id: asym.to_owned(),
            position: SourceResiduePosition::Author(sequence.to_owned()),
            insertion_code: None,
        };
        for (rows, links, breaks) in [
            // Each copy is linked, not only the first one listed.
            (
                [copy_b, copy_b_gal, copy_c, copy_c_gal].concat(),
                vec![("B", "B"), ("C", "C")],
                vec![],
            ),
            // B's unmodelled NAG leaves its GAL unlinked; C's NAG does not
            // stand in for it.
            (
                [copy_b_gal, copy_c, copy_c_gal].concat(),
                vec![("C", "C")],
                vec![],
            ),
            // A link beyond bonding distance is a reported break.
            (
                [copy_b, copy_b_gal, copy_c, copy_c_gal_far].concat(),
                vec![("B", "B")],
                vec![StructureInterpretIssue::ChainBreak {
                    residue: residue("C", "20"),
                    next: residue("C", "21"),
                }],
            ),
        ] {
            let text = two_branch_copies(&rows);
            let interpretation = interpret(&text);
            let report = interpretation.report();
            let asym_of = |atom| {
                report
                    .instances()
                    .iter()
                    .flat_map(|instance| instance.atoms())
                    .find(|candidate| candidate.atom() == atom)
                    .unwrap()
                    .asym_id()
                    .to_owned()
            };
            let mut actual = report
                .instances()
                .iter()
                .flat_map(|instance| instance.bonds())
                .filter(|bond| *bond.evidence() == BondEvidence::BranchLink)
                .map(|bond| bond.atoms().map(asym_of))
                .collect::<Vec<_>>();
            actual.sort();
            let expected = links
                .iter()
                .map(|(left, right)| [left.to_string(), right.to_string()])
                .collect::<Vec<_>>();
            assert_eq!(actual, expected, "{rows}");
            let actual_breaks = report
                .issues()
                .iter()
                .filter(|issue| matches!(issue, StructureInterpretIssue::ChainBreak { .. }))
                .cloned()
                .collect::<Vec<_>>();
            assert_eq!(actual_breaks, breaks, "{rows}");
        }
    }

    use crate::components::ComponentLibrary;
    use crate::macromolecular::{
        BondEvidence, ConnectivityOptions, SourceResidueId, SourceResiduePosition,
        StructureInterpretIssue, StructureInterpretOptions, StructureInterpretation,
    };

    /// A minimal block whose `rows` are `type name component asym seq x y z`.
    /// Coordinates of standard residues are the dictionary's ideal ones.
    fn block(rows: &str, extra: &str) -> String {
        records("HETATM", rows, extra)
    }

    /// Like [`block`]; `ATOM` records make the chain an inferred polymer.
    fn records(group: &str, rows: &str, extra: &str) -> String {
        let mut text = format!(
            "data_test\n{extra}\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
        );
        for (index, row) in rows
            .lines()
            .map(str::trim)
            .filter(|row| !row.is_empty())
            .enumerate()
        {
            text.push_str(&format!("{group} {} {row}\n", index + 1));
        }
        text
    }

    fn interpret(text: &str) -> StructureInterpretation {
        crate::mmcif::parse_str(text, Default::default())
            .unwrap()
            .interpret(Default::default())
            .unwrap()
    }

    /// Report issues other than the entity kinds these minimal blocks omit.
    fn issues(result: &StructureInterpretation) -> Vec<StructureInterpretIssue> {
        result
            .report()
            .issues()
            .iter()
            .filter(|issue| !matches!(issue, StructureInterpretIssue::EntityTypeInferred { .. }))
            .cloned()
            .collect()
    }

    fn line_of(text: &str, needle: &str) -> usize {
        text.lines().position(|line| line.contains(needle)).unwrap() + 1
    }

    /// Bonds by atom names, order, and evidence.
    fn bonds(result: &StructureInterpretation) -> Vec<(String, String, BondOrder, BondEvidence)> {
        let name = |atom| {
            result
                .report()
                .instances()
                .iter()
                .flat_map(|instance| instance.atoms())
                .find(|candidate| candidate.atom() == atom)
                .unwrap()
                .atom_name()
                .to_owned()
        };
        result
            .report()
            .instances()
            .iter()
            .flat_map(|instance| instance.bonds())
            .map(|bond| {
                let [a, b] = bond.atoms();
                let order = result
                    .topology()
                    .instance(a.instance())
                    .unwrap()
                    .molecule()
                    .bond_between(a.atom(), b.atom())
                    .map(|id| {
                        result
                            .topology()
                            .instance(a.instance())
                            .unwrap()
                            .molecule()
                            .bond(id)
                            .unwrap()
                            .order
                    })
                    .unwrap();
                (name(a), name(b), order, bond.evidence().clone())
            })
            .collect()
    }

    fn order_between(result: &StructureInterpretation, a: &str, b: &str) -> BondOrder {
        bonds(result)
            .into_iter()
            .find(|(left, right, _, _)| (left == a && right == b) || (left == b && right == a))
            .unwrap()
            .2
    }

    /// Nonzero formal charges by atom name.
    fn charges(result: &StructureInterpretation) -> Vec<(String, i8)> {
        result
            .report()
            .instances()
            .iter()
            .flat_map(|instance| instance.atoms())
            .filter_map(|atom| {
                let charge = result
                    .topology()
                    .instance(atom.atom().instance())
                    .unwrap()
                    .molecule()
                    .atom(atom.atom().atom())
                    .unwrap()
                    .formal_charge;
                (charge != 0).then(|| (atom.atom_name().to_owned(), charge))
            })
            .collect()
    }

    fn residue(asym: &str, sequence: i32) -> SourceResidueId {
        SourceResidueId {
            model_id: "1".to_owned(),
            asym_id: asym.to_owned(),
            position: SourceResiduePosition::Label(sequence),
            insertion_code: None,
        }
    }

    fn template(component: &str, variant: Option<&str>) -> BondEvidence {
        BondEvidence::Template {
            component_id: component.to_owned(),
            variant: variant.map(str::to_owned),
        }
    }

    const ASPARTATE: &str = "
        N N ASP A 1 -0.317 1.688 0.066
        C CA ASP A 1 -0.470 0.286 -0.344
        C C ASP A 1 -1.868 -0.180 -0.029
        O O ASP A 1 -2.534 0.415 0.786
        C CB ASP A 1 0.539 -0.580 0.413
        C CG ASP A 1 1.938 -0.195 0.004
        O OD1 ASP A 1 2.109 0.681 -0.810
        O OD2 ASP A 1 2.992 -0.826 0.543
        O OXT ASP A 1 -2.374 -1.256 -0.652
        H H ASP A 1 -0.928 2.289 -0.467
        H H2 ASP A 1 -0.478 1.795 1.056
        H HA ASP A 1 -0.292 0.199 -1.416
        H HB2 ASP A 1 0.419 -0.425 1.485
        H HB3 ASP A 1 0.367 -1.630 0.176
        H HD2 ASP A 1 3.869 -0.545 0.250
        H HXT ASP A 1 -3.275 -1.517 -0.416";

    const HISTIDINE: &str = "
        N N HIS A 1 -0.040 -1.210 0.053
        C CA HIS A 1 1.172 -1.709 0.652
        C C HIS A 1 1.083 -3.207 0.905
        O O HIS A 1 0.040 -3.770 1.222
        C CB HIS A 1 1.484 -0.975 1.962
        C CG HIS A 1 2.940 -1.060 2.353
        N ND1 HIS A 1 3.380 -2.075 3.129
        C CD2 HIS A 1 3.960 -0.251 2.046
        C CE1 HIS A 1 4.693 -1.908 3.317
        N NE2 HIS A 1 5.058 -0.801 2.662
        O OXT HIS A 1 2.247 -3.882 0.744
        H H HIS A 1 -0.102 -1.155 -0.950
        H H2 HIS A 1 -0.715 -0.741 0.634
        H HA HIS A 1 1.965 -1.558 -0.089
        H HB2 HIS A 1 1.215 0.087 1.879
        H HB3 HIS A 1 0.859 -1.368 2.775
        H HD1 HIS A 1 2.828 -2.838 3.511
        H HD2 HIS A 1 4.108 0.647 1.479
        H HE1 HIS A 1 5.340 -2.550 3.892
        H HE2 HIS A 1 6.002 -0.428 2.627
        H HXT HIS A 1 2.188 -4.848 0.901";

    fn without(rows: &str, atoms: &[&str]) -> String {
        rows.lines()
            .filter(|row| {
                !atoms
                    .iter()
                    .any(|atom| row.split_whitespace().nth(1) == Some(atom))
            })
            .collect::<Vec<_>>()
            .join("\n")
    }

    #[test]
    fn unstated_coordination_in_a_component_table_follows_template_then_elements() {
        // wwPDB entry files list `_chem_comp_bond` without the coordination
        // flag; the SF4 template flags its Fe-S bonds, so they are dative as
        // a PDB reading of the same structure makes them.
        let table = |flag_column: bool, rows: &str| {
            let mut text = String::from(
                "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n",
            );
            if flag_column {
                text.push_str("_chem_comp_bond.pdbx_metal_coordination_flag\n");
            }
            text.push_str(rows);
            text
        };
        let sf4 = "FE FE1 SF4 A 1 40.988 0.074 22.637
             S S2 SF4 A 1 43.144 0.781 22.394";
        let unstated = interpret(&block(sf4, &table(false, "SF4 FE1 S2 sing\n")));
        assert_eq!(
            bonds(&unstated)
                .iter()
                .map(|(_, _, order, evidence)| (*order, evidence.clone()))
                .collect::<Vec<_>>(),
            [(BondOrder::Dative, BondEvidence::ComponentTable)]
        );
        // A stated flag wins over the template.
        let stated = interpret(&block(sf4, &table(true, "SF4 FE1 S2 sing N\n")));
        assert_eq!(bonds(&stated)[0].2, BondOrder::Single);
        // Without a template, a metal and a non-metal are coordinated, as for
        // PDB `LINK` and `CONECT`; other pairs keep their stated order.
        let unknown = interpret(&block(
            "FE FE XFC A 1 0.000 0.000 0.000
             CL CL1 XFC A 1 2.200 0.000 0.000
             C C1 XFC A 1 -1.900 0.000 0.000
             O O1 XFC A 1 -3.050 0.000 0.000",
            &table(false, "XFC FE CL1 sing\nXFC C1 O1 doub\n"),
        ));
        let orders = bonds(&unknown)
            .into_iter()
            .map(|(a, b, order, _)| (a, b, order))
            .collect::<Vec<_>>();
        assert!(orders.contains(&("FE".to_owned(), "CL1".to_owned(), BondOrder::Dative)));
        assert!(orders.contains(&("C1".to_owned(), "O1".to_owned(), BondOrder::Double)));
    }

    #[test]
    fn templates_supply_bonds_and_metal_coordination_is_dative() {
        // The SF4 template flags every Fe-S bond as metal coordination;
        // its sulfides carry -2 and its iron charges are unknown.
        let cluster = interpret(&block(
            "FE FE1 SF4 A 1 40.988 0.074 22.637
             FE FE2 SF4 A 1 39.644 2.319 22.749
             FE FE3 SF4 A 1 42.036 2.170 23.811
             FE FE4 SF4 A 1 41.755 2.132 21.208
             S S1 SF4 A 1 41.267 3.911 22.551
             S S2 SF4 A 1 43.144 0.781 22.394
             S S3 SF4 A 1 39.806 0.980 20.921
             S S4 SF4 A 1 40.200 1.034 24.551",
            "",
        ));
        assert_eq!(cluster.topology().instance_count(), 1);
        let cluster_bonds = bonds(&cluster);
        assert_eq!(cluster_bonds.len(), 12);
        assert!(cluster_bonds
            .iter()
            .all(|(_, _, order, evidence)| *order == BondOrder::Dative
                && *evidence == template("SF4", None)));
        assert_eq!(
            charges(&cluster),
            [("S1", -2), ("S2", -2), ("S3", -2), ("S4", -2)]
                .map(|(name, charge)| (name.to_owned(), charge))
        );

        // A metal-coordination link joins a zinc ion and a cysteine into one
        // molecule; the ion's charge comes from its template.
        let zinc = interpret(&block(
            "N N CYS A 1 1.585 0.483 -0.081
             C CA CYS A 1 0.141 0.450 0.186
             C C CYS A 1 -0.095 0.006 1.606
             O O CYS A 1 0.685 -0.742 2.143
             C CB CYS A 1 -0.533 -0.530 -0.774
             S SG CYS A 1 -0.247 0.004 -2.484
             ZN ZN ZN B . -0.247 0.004 -4.784",
            "loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_label_seq_id\nmetalc1 metalc A SG 1 B ZN .\n",
        ));
        assert_eq!(zinc.topology().instance_count(), 1);
        let cysteine = template("CYS", None);
        let as_strings = |rows: &[(&str, &str, BondOrder, BondEvidence)]| {
            rows.iter()
                .map(|(a, b, order, evidence)| {
                    ((*a).to_owned(), (*b).to_owned(), *order, evidence.clone())
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(
            bonds(&zinc),
            as_strings(&[
                (
                    "SG",
                    "ZN",
                    BondOrder::Dative,
                    BondEvidence::Connection {
                        connection_id: Some("metalc1".to_owned()),
                        connection_type: "metalc".to_owned(),
                    },
                ),
                ("N", "CA", BondOrder::Single, cysteine.clone()),
                ("CA", "C", BondOrder::Single, cysteine.clone()),
                ("CA", "CB", BondOrder::Single, cysteine.clone()),
                ("C", "O", BondOrder::Double, cysteine.clone()),
                ("CB", "SG", BondOrder::Single, cysteine),
            ])
        );
        assert_eq!(charges(&zinc), [("ZN".to_owned(), 2)]);
    }

    #[test]
    fn model_hydrogens_set_heteroatom_charges_and_select_tautomers() {
        // Without model hydrogens the template's charges apply: the
        // dictionary's histidine is the ND1-protonated imidazolium.
        let heavy = interpret(&block(
            &without(
                HISTIDINE,
                &[
                    "H", "H2", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "HXT",
                ],
            ),
            "",
        ));
        assert_eq!(charges(&heavy), [("ND1".to_owned(), 1)]);
        assert!(issues(&heavy).is_empty());

        // An extra amine hydrogen named outside the template, or outside the
        // block's own bond table, attaches by distance; with OXT unprotonated
        // the aspartate is a zwitterion.
        let table = "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n\
            _chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n\
            ASP N CA sing\nASP CA C sing\nASP C O doub\nASP C OXT sing\nASP CA CB sing\n\
            ASP CB CG sing\nASP CG OD1 doub\nASP CG OD2 sing\nASP N H sing\nASP N H2 sing\n\
            ASP CA HA sing\nASP CB HB2 sing\nASP CB HB3 sing\nASP OD2 HD2 sing\n";
        for extra in ["", table] {
            let zwitterion = interpret(&block(
                &format!(
                    "{}\nH H3 ASP A 1 0.633 1.938 -0.134",
                    without(ASPARTATE, &["HXT"])
                ),
                extra,
            ));
            assert_eq!(zwitterion.topology().instance_count(), 1, "{extra}");
            assert_eq!(
                charges(&zwitterion),
                [("N".to_owned(), 1), ("OXT".to_owned(), -1)],
                "{extra}"
            );
            assert_eq!(
                issues(&zwitterion),
                [
                    StructureInterpretIssue::HydrogenAttachedByDistance {
                        residue: residue("A", 1),
                        atom_name: "H3".to_owned(),
                        heavy_atom_name: "N".to_owned(),
                    },
                    StructureInterpretIssue::ProtonationFromHydrogens {
                        residue: residue("A", 1),
                        component_id: "ASP".to_owned(),
                        variant: None,
                        charged_atoms: vec![("N".to_owned(), 1), ("OXT".to_owned(), -1)],
                    },
                ],
                "{extra}"
            );
        }

        // A deprotonated side chain carries the charge on its acid oxygen.
        let aspartate = interpret(&block(&without(ASPARTATE, &["HD2"]), ""));
        assert_eq!(charges(&aspartate), [("OD2".to_owned(), -1)]);

        // A ring nitrogen takes or loses one proton whatever its bond orders:
        // the cytosine base (ideal coordinates of the dictionary's DC)
        // protonated at N3, as in i-motif DNA, and the imidazolate histidine.
        let cytosine = "
            N N1 DC A 1 -0.267 -0.584 2.528
            C C2 DC A 1 0.270 0.648 2.563
            O O2 DC A 1 0.052 1.424 1.647
            N N3 DC A 1 1.037 1.035 3.581
            C C4 DC A 1 1.291 0.212 4.589
            N N4 DC A 1 2.085 0.622 5.635
            C C5 DC A 1 0.746 -1.088 4.580
            C C6 DC A 1 -0.035 -1.465 3.541
            H H41 DC A 1 2.461 1.516 5.636
            H H42 DC A 1 2.265 0.019 6.373
            H H5 DC A 1 0.943 -1.771 5.394
            H H6 DC A 1 -0.467 -2.454 3.514
            H H3 DC A 1 1.431 2.041 3.589";
        let imidazolate = without(HISTIDINE, &["HD1", "HE2"]);
        for (rows, component, charged, attached) in [
            (cytosine, "DC", vec![("N3".to_owned(), 1)], Some("N3")),
            (
                imidazolate.as_str(),
                "HIS",
                vec![("ND1".to_owned(), 0), ("NE2".to_owned(), -1)],
                None,
            ),
        ] {
            let result = interpret(&block(rows, ""));
            let nonzero = charged
                .iter()
                .filter(|(_, charge)| *charge != 0)
                .cloned()
                .collect::<Vec<_>>();
            assert_eq!(charges(&result), nonzero, "{component}");
            let mut expected = attached
                .map(
                    |heavy| StructureInterpretIssue::HydrogenAttachedByDistance {
                        residue: residue("A", 1),
                        atom_name: "H3".to_owned(),
                        heavy_atom_name: heavy.to_owned(),
                    },
                )
                .into_iter()
                .collect::<Vec<_>>();
            expected.push(StructureInterpretIssue::ProtonationFromHydrogens {
                residue: residue("A", 1),
                component_id: component.to_owned(),
                variant: None,
                charged_atoms: charged,
            });
            assert_eq!(issues(&result), expected, "{component}");
        }

        // Missing ring hydrogens select a neutral tautomer, whose ring bonds
        // follow the hydrogen.
        for (removed, variant, nd1_ce1, ce1_ne2) in [
            (
                "HD1",
                "NE2-protonated",
                BondOrder::Double,
                BondOrder::Single,
            ),
            (
                "HE2",
                "ND1-protonated",
                BondOrder::Single,
                BondOrder::Double,
            ),
        ] {
            let tautomer = interpret(&block(&without(HISTIDINE, &[removed]), ""));
            assert_eq!(charges(&tautomer), [], "{variant}");
            assert_eq!(
                issues(&tautomer),
                [StructureInterpretIssue::ProtonationFromHydrogens {
                    residue: residue("A", 1),
                    component_id: "HIS".to_owned(),
                    variant: Some(variant.to_owned()),
                    charged_atoms: Vec::new(),
                }]
            );
            assert_eq!(order_between(&tautomer, "ND1", "CE1"), nd1_ce1, "{variant}");
            assert_eq!(order_between(&tautomer, "CE1", "NE2"), ce1_ne2, "{variant}");
        }
    }

    const UNKNOWN_LIGAND: &str = "
        C C1 XYZ A 1 0.0 0.0 0.0
        C C2 XYZ A 1 1.5 0.0 0.0
        O O1 XYZ A 1 2.0 1.2 0.0";

    #[test]
    fn residues_without_templates_infer_bonds_unless_disabled() {
        let text = block(UNKNOWN_LIGAND, "");
        let inferred = interpret(&text);
        assert_eq!(
            bonds(&inferred),
            [
                ("C1", "C2", BondOrder::Single, BondEvidence::Distance),
                ("C2", "O1", BondOrder::Single, BondEvidence::Distance),
            ]
            .map(|(a, b, order, evidence)| (
                a.to_owned(),
                b.to_owned(),
                order,
                evidence
            ))
        );
        assert_eq!(
            issues(&inferred),
            [StructureInterpretIssue::ResidueBondsInferred {
                residue: residue("A", 1),
                component_id: "XYZ".to_owned(),
                bonds: 2,
            }]
        );

        let strict = StructureInterpretOptions::default().with_connectivity(
            ConnectivityOptions::default().with_infer_bonds_from_distance(false),
        );
        let error = crate::mmcif::parse_str(&text, Default::default())
            .unwrap()
            .interpret(strict.clone())
            .unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::Unsupported);
        assert_eq!(
            error.location().map(|location| location.line()),
            Some(line_of(&text, " C1 XYZ"))
        );

        // A caller template takes precedence and supplies the bond order.
        let library = ComponentLibrary::from_mmcif(
            &crate::mmcif::parse_str(
                "data_XYZ\n_chem_comp.id XYZ\n_chem_comp.type NON-POLYMER\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\nXYZ C1 C 0\nXYZ C2 C 0\nXYZ O1 O 0\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nXYZ C1 C2 SING\nXYZ C2 O1 DOUB\n",
                Default::default(),
            )
            .unwrap(),
        )
        .unwrap();
        let with_library = strict.clone().with_connectivity(
            strict
                .connectivity()
                .clone()
                .with_components(Some(std::sync::Arc::new(library))),
        );
        let templated = crate::mmcif::parse_str(&text, Default::default())
            .unwrap()
            .interpret(with_library)
            .unwrap();
        assert_eq!(
            bonds(&templated),
            [
                ("C1", "C2", BondOrder::Single, template("XYZ", None)),
                ("C2", "O1", BondOrder::Double, template("XYZ", None)),
            ]
            .map(|(a, b, order, evidence)| (
                a.to_owned(),
                b.to_owned(),
                order,
                evidence
            ))
        );
        assert!(issues(&templated).is_empty());
    }

    const DIGLYCINE: &str = "
        N N GLY A 1 0.0 0.0 0.0
        C CA GLY A 1 1.45 0.0 0.0
        C C GLY A 1 2.0 1.4 0.0
        O O GLY A 1 1.3 2.4 0.0
        N N GLY A 2 3.33 1.4 0.0
        C CA GLY A 2 4.0 2.6 0.0
        C C GLY A 2 5.5 2.6 0.0
        O O GLY A 2 6.1 3.6 0.0";

    #[test]
    fn polymer_links_require_bonding_distance() {
        // Without _entity_poly, both residues' templates identify a peptide.
        let linked = interpret(&records("ATOM", DIGLYCINE, ""));
        assert_eq!(linked.topology().instance_count(), 1);
        assert_eq!(
            bonds(&linked)
                .iter()
                .filter(|(_, _, _, evidence)| *evidence == BondEvidence::PolymerLink)
                .map(|(a, b, _, _)| (a.clone(), b.clone()))
                .collect::<Vec<_>>(),
            [("C".to_owned(), "N".to_owned())]
        );
        assert!(issues(&linked).is_empty());

        // Residue 2 moved 3 A away: adjacent in sequence, but not bonded.
        let shifted = DIGLYCINE
            .lines()
            .map(|row| {
                let fields = row.split_whitespace().collect::<Vec<_>>();
                if fields.get(4) == Some(&"2") {
                    let x = fields[5].parse::<f64>().unwrap() + 3.0;
                    format!(
                        "{} {} {} {} {} {x} {} {}",
                        fields[0], fields[1], fields[2], fields[3], fields[4], fields[6], fields[7]
                    )
                } else {
                    row.to_owned()
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        let broken = interpret(&records("ATOM", &shifted, ""));
        assert_eq!(broken.topology().instance_count(), 2);
        assert_eq!(broken.topology().bond_count(), 6);
        assert_eq!(
            issues(&broken),
            [StructureInterpretIssue::ChainBreak {
                residue: residue("A", 1),
                next: residue("A", 2),
            }]
        );
    }

    #[test]
    fn untemplated_polymer_residues_link_by_their_neighbours_linkage_or_atom_names() {
        let inferred = |sequence, component: &str| StructureInterpretIssue::ResidueBondsInferred {
            residue: residue("A", sequence),
            component_id: component.to_owned(),
            bonds: 3,
        };
        // Without _entity_poly, a residue without a template defers to its
        // templated neighbour, and two such residues link by atom names.
        for (first, second, expected) in [
            ("GLY", "XGL", vec![inferred(2, "XGL")]),
            ("XGL", "GLY", vec![inferred(1, "XGL")]),
            ("XGL", "YGL", vec![inferred(1, "XGL"), inferred(2, "YGL")]),
        ] {
            let rows = DIGLYCINE
                .replace("GLY A 1", &format!("{first} A 1"))
                .replace("GLY A 2", &format!("{second} A 2"));
            let linked = interpret(&records("ATOM", &rows, ""));
            let polymer_links = bonds(&linked)
                .iter()
                .filter(|(_, _, _, evidence)| *evidence == BondEvidence::PolymerLink)
                .map(|(a, b, _, _)| (a.clone(), b.clone()))
                .collect::<Vec<_>>();
            assert_eq!(polymer_links, [("C".to_owned(), "N".to_owned())]);
            assert_eq!(linked.topology().instance_count(), 1);
            assert_eq!(issues(&linked), expected);
        }
    }

    #[test]
    fn template_mismatches_fail_and_split_residues_are_reported() {
        let unknown_atom = block(
            "N N GLY A 1 0.0 0.0 0.0
             C CA GLY A 1 1.45 0.0 0.0
             C CX GLY A 1 2.0 1.4 0.0",
            "",
        );
        let error = crate::mmcif::parse_str(&unknown_atom, Default::default())
            .unwrap()
            .interpret(Default::default())
            .unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
        assert_eq!(
            error.location().map(|location| location.line()),
            Some(line_of(&unknown_atom, " CX GLY"))
        );
        assert_eq!(error.message(), "atom `CX` is not in component `GLY`");

        // Without its CA, a glycine's amine is cut off from its carbonyl.
        let split = interpret(&block(
            "N N GLY A 1 0.0 0.0 0.0
             C C GLY A 1 2.0 1.4 0.0
             O O GLY A 1 1.3 2.4 0.0",
            "",
        ));
        assert_eq!(split.topology().instance_count(), 2);
        assert_eq!(
            issues(&split),
            [StructureInterpretIssue::ResidueFragmented {
                residue: residue("A", 1),
                fragments: 2,
            }]
        );
    }

    #[test]
    fn stated_internal_bonds_take_precedence_over_templates() {
        // Files that enumerate a residue's bonds, as Kekule's writer does,
        // are taken as stated; the template does not add the others.
        let text = block(
            &DIGLYCINE.lines().take(5).collect::<Vec<_>>().join("\n"),
            "loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_label_seq_id\nc1 covale A N 1 A CA 1\n",
        );
        let stated = interpret(&text);
        assert_eq!(
            bonds(&stated),
            [(
                "N".to_owned(),
                "CA".to_owned(),
                BondOrder::Single,
                BondEvidence::Connection {
                    connection_id: Some("c1".to_owned()),
                    connection_type: "covale".to_owned(),
                },
            )]
        );
    }
}
