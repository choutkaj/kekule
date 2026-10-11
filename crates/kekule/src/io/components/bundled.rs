//! Templates bundled with Kekule. Their provenance, selection, and checksums
//! are recorded in `tools/components`.

use super::super::{parse_mmcif_str, MmcifParseOptions};
use super::{Component, ComponentBondOrder, ComponentLibrary};

macro_rules! bundled_components {
    [$($id:literal),* $(,)?] => {
        &[$(($id, include_str!(concat!("../../../data/components/", $id, ".cif")))),*]
    };
}

/// Unmodified Chemical Component Dictionary entries, grouped as in
/// `tools/components/fetch.py`.
#[rustfmt::skip]
const ENTRIES: &[(&str, &str)] = bundled_components![
    // Amino acids.
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "SEC", "PYL", "MSE",
    // Nucleotides.
    "A", "C", "G", "U", "DA", "DC", "DG", "DT",
    // Water.
    "HOH", "DOD",
    // Ions.
    "ZN", "MG", "CA", "NA", "K", "CL", "MN", "FE", "FE2", "CU", "CU1",
    "NI", "CO", "CD", "HG", "IOD", "BR",
    // Cofactors.
    "HEM", "HEC", "FAD", "FMN", "NAD", "NAI", "NAP", "NDP", "ATP", "ADP",
    "AMP", "GTP", "GDP", "SAM", "SAH", "COA", "ACO", "PLP", "TPP", "BTN",
    "SF4", "FES", "F3S", "CLA",
    // Crystallization additives.
    "SO4", "PO4", "GOL", "EDO", "ACT", "PEG", "DMS", "MPD",
];

pub(super) fn library() -> ComponentLibrary {
    let mut library = ComponentLibrary::new();
    for (id, text) in ENTRIES {
        let document = parse_mmcif_str(text, MmcifParseOptions::default())
            .unwrap_or_else(|error| panic!("bundled component {id} parses: {error}"));
        let component = Component::from_mmcif_block(&document.blocks()[0])
            .unwrap_or_else(|error| panic!("bundled component {id} is valid: {error}"));
        library.insert(component);
    }
    let histidine = library.get("HIS").expect("HIS is bundled").clone();
    for variant in histidine_tautomers(&histidine) {
        library
            .insert_variant(variant)
            .expect("histidine variants name HIS");
    }
    library
}

/// The neutral histidine tautomers. The dictionary defines only the
/// imidazolium form, with both ring hydrogens and ND1 charged.
fn histidine_tautomers(histidine: &Component) -> [Component; 2] {
    // NE2 carries the ring hydrogen: drop HD1; the ND1=CE1 bond stays double.
    let mut ne2 = histidine.without_atom("HD1");
    ne2.set_formal_charge("ND1", 0);
    ne2.variant = Some("NE2-protonated".to_owned());
    // ND1 carries the ring hydrogen: drop HE2 and move the double bond to CE1=NE2.
    let mut nd1 = histidine.without_atom("HE2");
    nd1.set_formal_charge("ND1", 0);
    nd1.set_bond_order("ND1", "CE1", ComponentBondOrder::Single);
    nd1.set_bond_order("CE1", "NE2", ComponentBondOrder::Double);
    nd1.variant = Some("ND1-protonated".to_owned());
    [ne2, nd1]
}

impl Component {
    pub(super) fn without_atom(&self, name: &str) -> Self {
        let (removed, _) = self.atom(name).expect("template atom exists");
        let mut component = self.clone();
        component.atoms.remove(removed);
        component
            .bonds
            .retain(|bond| !bond.atoms.contains(&removed));
        for bond in &mut component.bonds {
            for atom in &mut bond.atoms {
                if *atom > removed {
                    *atom -= 1;
                }
            }
        }
        component
    }

    pub(super) fn set_formal_charge(&mut self, name: &str, charge: i8) {
        let (index, _) = self.atom(name).expect("template atom exists");
        self.atoms[index].formal_charge = Some(charge);
    }

    /// Renames atoms all at once, so that `[("HB2", "HB1"), ("HB3", "HB2")]`
    /// shifts names rather than chaining them. Alternate names are cleared.
    pub(super) fn rename_atoms(&mut self, renamed: &[(&str, &str)]) {
        let indices = renamed
            .iter()
            .map(|&(name, _)| self.atom(name).expect("template atom exists").0)
            .collect::<Vec<_>>();
        for (index, &(_, name)) in indices.into_iter().zip(renamed) {
            self.atoms[index].name = name.to_owned();
            self.atoms[index].alternate_name = None;
        }
    }

    fn set_bond_order(&mut self, a: &str, b: &str, order: ComponentBondOrder) {
        let (a, _) = self.atom(a).expect("template atom exists");
        let (b, _) = self.atom(b).expect("template atom exists");
        let bond = self
            .bonds
            .iter_mut()
            .find(|bond| bond.atoms == [a, b] || bond.atoms == [b, a])
            .expect("template bond exists");
        bond.order = order;
    }
}
