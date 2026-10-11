//! Templates for the residue names of molecular-mechanics force fields,
//! derived from the bundled dictionary entries.
//!
//! Force fields name protonation states, termini, water, and ions with
//! residue and atom names of their own, such as AMBER `HIE` for the
//! NE2-protonated histidine or CHARMM `HSD` for the ND1-protonated one. Each
//! template here is a bundled entry renamed, with hydrogens added or removed
//! and formal charges set as the force field defines the residue. Names follow
//! the residue definitions of the AMBER ff14SB, OL15, and OL3 and the CHARMM36
//! libraries.

use super::{Component, ComponentAtom, ComponentBond, ComponentBondOrder, ComponentLibrary};
use crate::core::Element;

/// Amino acids AMBER defines N- and C-terminal residues for, as `NALA` and
/// `CALA`.
const AMBER_TERMINAL: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU", "GLY", "HID", "HIE", "HIP", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];

pub(super) fn amber() -> ComponentLibrary {
    let bundled = ComponentLibrary::bundled();
    let base = |id: &str| bundled.get(id).expect("bundled component exists").clone();
    let histidine = |variant: &str| histidine(bundled, variant).clone();
    let mut residues = vec![
        histidine("ND1-protonated").named("HID", "histidine, ND1-protonated"),
        histidine("NE2-protonated").named("HIE", "histidine, NE2-protonated"),
        base("HIS").named("HIP", "histidine, doubly protonated"),
        base("CYS")
            .without_atom("HG")
            .named("CYX", "cysteine in a disulfide"),
        base("CYS")
            .without_atom("HG")
            .charged("SG", -1)
            .named("CYM", "cysteine, deprotonated"),
        base("ASP")
            .without_atom("HD2")
            .charged("OD2", -1)
            .named("ASP", "aspartate"),
        base("ASP").named("ASH", "aspartic acid"),
        base("GLU")
            .without_atom("HE2")
            .charged("OE2", -1)
            .named("GLU", "glutamate"),
        base("GLU").named("GLH", "glutamic acid"),
        base("LYS")
            .without_atom("HZ1")
            .charged("NZ", 0)
            .named("LYN", "lysine, neutral"),
        // The acetyl cap is acetate without its leaving oxygen; ff14SB names
        // its methyl hydrogens HH31-HH33, later libraries H1-H3.
        base("ACT")
            .without_atom("OXT")
            .alternates(&[("H1", "HH31"), ("H2", "HH32"), ("H3", "HH33")])
            .named("ACE", "acetyl cap"),
        base("HOH").named("WAT", "water"),
        // Ion names of older libraries; current ones use dictionary names.
        base("NA")
            .renamed(&[("NA", "Na+")])
            .named("Na+", "sodium ion"),
        base("K")
            .renamed(&[("K", "K+")])
            .named("K+", "potassium ion"),
        base("CL")
            .renamed(&[("CL", "Cl-")])
            .named("Cl-", "chloride ion"),
    ];
    // Nucleotides keep the atom names of their residues. A linking phosphate
    // is an anion; 5'-terminal and free nucleotides have none, and atoms a
    // residue lacks are simply absent.
    for (nucleotide, name) in [
        ("DA", "deoxyadenosine"),
        ("DC", "deoxycytidine"),
        ("DG", "deoxyguanosine"),
        ("DT", "thymidine"),
        ("A", "adenosine"),
        ("C", "cytidine"),
        ("G", "guanosine"),
        ("U", "uridine"),
    ] {
        let phosphate = || base(nucleotide).without_atom("HOP2").charged("OP2", -1);
        residues.push(phosphate().named(nucleotide, &format!("{name} monophosphate")));
        residues.push(phosphate().named(
            &format!("{nucleotide}3"),
            &format!("3'-terminal {name} monophosphate"),
        ));
        residues.push(
            base(nucleotide).named(&format!("{nucleotide}5"), &format!("5'-terminal {name}")),
        );
        residues.push(base(nucleotide).named(&format!("{nucleotide}N"), &format!("free {name}")));
    }
    let mut library = ComponentLibrary::new();
    for residue in residues {
        library.insert(residue);
    }
    // N-terminal residues are ammonium ions with hydrogens H1-H3 (H2-H3 on
    // proline); C-terminal ones are carboxylates.
    for &id in AMBER_TERMINAL {
        let residue = library.get(id).cloned().unwrap_or_else(|| base(id));
        let name = describe(&residue);
        let ammonium = if id == "PRO" {
            residue.clone().renamed(&[("H", "H2")])
        } else {
            residue.clone().alternates(&[("H", "H1")])
        };
        library.insert(
            ammonium
                .with_hydrogen("H3", "N")
                .charged("N", 1)
                .named(&format!("N{id}"), &format!("N-terminal {name}")),
        );
        library.insert(
            residue
                .without_atom("HXT")
                .charged("OXT", -1)
                .named(&format!("C{id}"), &format!("C-terminal {name}")),
        );
    }
    library.retitle("AMBER")
}

/// Hydrogen renames from dictionary to CHARMM names: methylene hydrogens are
/// numbered 1-2 rather than 2-3.
#[rustfmt::skip]
const CHARMM_RENAMES: &[(&str, &[(&str, &str)])] = &[
    ("ALA", &[]),
    ("ARG", &[("HB2", "HB1"), ("HB3", "HB2"), ("HG2", "HG1"), ("HG3", "HG2"), ("HD2", "HD1"), ("HD3", "HD2")]),
    ("ASN", &[("HB2", "HB1"), ("HB3", "HB2")]),
    ("ASP", &[("HB2", "HB1"), ("HB3", "HB2")]),
    ("CYS", &[("HB2", "HB1"), ("HB3", "HB2"), ("HG", "HG1")]),
    ("GLN", &[("HB2", "HB1"), ("HB3", "HB2"), ("HG2", "HG1"), ("HG3", "HG2")]),
    ("GLU", &[("HB2", "HB1"), ("HB3", "HB2"), ("HG2", "HG1"), ("HG3", "HG2")]),
    ("GLY", &[("HA2", "HA1"), ("HA3", "HA2")]),
    ("HIS", &[("HB2", "HB1"), ("HB3", "HB2")]),
    ("ILE", &[("HG12", "HG11"), ("HG13", "HG12"), ("CD1", "CD"), ("HD11", "HD1"), ("HD12", "HD2"), ("HD13", "HD3")]),
    ("LEU", &[("HB2", "HB1"), ("HB3", "HB2")]),
    ("LYS", &[("HB2", "HB1"), ("HB3", "HB2"), ("HG2", "HG1"), ("HG3", "HG2"), ("HD2", "HD1"), ("HD3", "HD2"), ("HE2", "HE1"), ("HE3", "HE2")]),
    ("MET", &[("HB2", "HB1"), ("HB3", "HB2"), ("HG2", "HG1"), ("HG3", "HG2")]),
    ("PHE", &[("HB2", "HB1"), ("HB3", "HB2")]),
    ("PRO", &[("HB2", "HB1"), ("HB3", "HB2"), ("HG2", "HG1"), ("HG3", "HG2"), ("HD2", "HD1"), ("HD3", "HD2")]),
    ("SER", &[("HB2", "HB1"), ("HB3", "HB2"), ("HG", "HG1")]),
    ("THR", &[]),
    ("TRP", &[("HB2", "HB1"), ("HB3", "HB2")]),
    ("TYR", &[("HB2", "HB1"), ("HB3", "HB2")]),
    ("VAL", &[]),
];

pub(super) fn charmm() -> ComponentLibrary {
    let bundled = ComponentLibrary::bundled();
    let base = |id: &str| bundled.get(id).expect("bundled component exists");
    let amino_acid = |id: &str, from: &Component| {
        let (_, renames) = CHARMM_RENAMES
            .iter()
            .find(|(residue, _)| *residue == id)
            .expect("CHARMM renames exist");
        // Backbone: HN for the amide hydrogen, and the terminal patches'
        // OT1 and OT2 for the carboxyl oxygens. Proline has no amide
        // hydrogen; its N-terminal patch names the ammonium hydrogens HN1 and
        // HN2. Terminal hydrogens HT1-HT3 name different atoms in different
        // patches, so they attach by distance.
        let amide = if id == "PRO" { "HN1" } else { "HN" };
        from.clone()
            .renamed(&[*renames, &[("H", amide), ("OXT", "OT2")]].concat())
            .alternates(&[("O", "OT1")])
    };
    let mut library = ComponentLibrary::new();
    for &(id, _) in CHARMM_RENAMES {
        if id != "HIS" {
            let residue = amino_acid(id, base(id));
            let name = describe(&residue);
            library.insert(residue.named(id, &name));
        }
    }
    let residues = [
        amino_acid("HIS", histidine(bundled, "ND1-protonated"))
            .named("HSD", "histidine, ND1-protonated"),
        amino_acid("HIS", histidine(bundled, "NE2-protonated"))
            .named("HSE", "histidine, NE2-protonated"),
        amino_acid("HIS", base("HIS")).named("HSP", "histidine, doubly protonated"),
        amino_acid("ASP", base("ASP"))
            .without_atom("HD2")
            .charged("OD2", -1)
            .named("ASP", "aspartate"),
        amino_acid("ASP", base("ASP")).named("ASPP", "aspartic acid"),
        amino_acid("GLU", base("GLU"))
            .without_atom("HE2")
            .charged("OE2", -1)
            .named("GLU", "glutamate"),
        amino_acid("GLU", base("GLU")).named("GLUP", "glutamic acid"),
        amino_acid("LYS", base("LYS"))
            .without_atom("HZ3")
            .charged("NZ", 0)
            .named("LSN", "lysine, neutral"),
        amino_acid("CYS", base("CYS"))
            .without_atom("HG1")
            .charged("SG", -1)
            .named("CYM", "cysteine, deprotonated"),
        base("HOH")
            .clone()
            .renamed(&[("O", "OH2")])
            .named("TIP3", "water"),
        base("NA")
            .clone()
            .renamed(&[("NA", "SOD")])
            .named("SOD", "sodium ion"),
        base("K")
            .clone()
            .renamed(&[("K", "POT")])
            .named("POT", "potassium ion"),
        base("CA")
            .clone()
            .renamed(&[("CA", "CAL")])
            .named("CAL", "calcium ion"),
        base("CL")
            .clone()
            .renamed(&[("CL", "CLA")])
            .named("CLA", "chloride ion"),
        base("ZN").clone().named("ZN2", "zinc ion"),
        base("CD").clone().named("CD2", "cadmium ion"),
    ];
    for residue in residues {
        library.insert(residue);
    }
    library.retitle("CHARMM")
}

fn histidine<'a>(library: &'a ComponentLibrary, variant: &str) -> &'a Component {
    library
        .variants("HIS")
        .find(|candidate| candidate.variant() == Some(variant))
        .expect("bundled histidine variant exists")
}

/// A lower-case description: dictionary names are upper case.
fn describe(component: &Component) -> String {
    if component
        .name
        .chars()
        .any(|character| character.is_ascii_lowercase())
    {
        component.name.clone()
    } else {
        component.name.to_ascii_lowercase()
    }
}

impl Component {
    /// The template under a force-field residue name, no longer a variant.
    fn named(mut self, id: &str, description: &str) -> Self {
        self.id = id.to_owned();
        self.name = description.to_owned();
        self.variant = None;
        self
    }

    fn charged(mut self, atom: &str, charge: i8) -> Self {
        self.set_formal_charge(atom, charge);
        self
    }

    fn renamed(mut self, renamed: &[(&str, &str)]) -> Self {
        self.rename_atoms(renamed);
        self
    }

    /// Gives atoms a further name, as `(name, alternate)`.
    fn alternates(mut self, alternates: &[(&str, &str)]) -> Self {
        for &(atom, alternate) in alternates {
            let (index, _) = self.atom(atom).expect("template atom exists");
            self.atoms[index].alternate_name = Some(alternate.to_owned());
        }
        self
    }

    /// Adds a hydrogen singly bonded to `heavy`.
    fn with_hydrogen(mut self, name: &str, heavy: &str) -> Self {
        let (heavy, _) = self.atom(heavy).expect("template atom exists");
        self.atoms.push(ComponentAtom {
            name: name.to_owned(),
            alternate_name: None,
            element: Element::from_atomic_number(1).expect("hydrogen"),
            isotope: None,
            formal_charge: Some(0),
            aromatic: false,
            leaving: false,
        });
        self.bonds.push(ComponentBond {
            atoms: [heavy, self.atoms.len() - 1],
            order: ComponentBondOrder::Single,
            aromatic: false,
            metal_coordination: false,
        });
        self
    }
}

impl ComponentLibrary {
    /// Names each template after the force field that defines its residue.
    fn retitle(self, force_field: &str) -> Self {
        let mut library = Self::new();
        for component in self.components() {
            let mut component = component.clone();
            component.name = format!("{} ({force_field} {})", component.name, component.id);
            library.insert(component);
        }
        library
    }
}
