//! Structural mmCIF parsing, interpretation, and writing.
//!
//! Data blocks are independent interpretation scopes. One selected coordinate
//! model naturally produces a [`crate::structure::Structure`], while compatible
//! coordinate models from one block may produce a
//! [`crate::structure::StructureSet`]. Ordinary writing derives mmCIF entity kinds
//! from canonical topology classification; source reports and explicit
//! classifications remain available for faithful round trips and expert
//! overrides.
//!
//! Interpretation keeps selected `_atom_site` rows in source order as dense atom
//! order, and writers emit `_atom_site` rows in dense atom order, so coordinate
//! files in the same order address the same atoms after a round trip.
//!
//! ```
//! use kekule::mmcif;
//!
//! let text = "\
//! data_water
//! loop_
//! _atom_site.group_PDB
//! _atom_site.id
//! _atom_site.type_symbol
//! _atom_site.label_atom_id
//! _atom_site.label_comp_id
//! _atom_site.label_asym_id
//! _atom_site.label_entity_id
//! _atom_site.label_seq_id
//! _atom_site.Cartn_x
//! _atom_site.Cartn_y
//! _atom_site.Cartn_z
//! _atom_site.auth_seq_id
//! _atom_site.auth_asym_id
//! _atom_site.pdbx_PDB_model_num
//! HETATM 1 O O   HOH A 1 . 0.000 0.000 0.000 1 A 1
//! HETATM 2 H H1  HOH A 1 . 0.960 0.000 0.000 1 A 1
//! HETATM 3 H H2  HOH A 1 . -0.240 0.930 0.000 1 A 1
//! ";
//! let document = mmcif::parse_str(text, Default::default())?;
//! let interpretation = document.interpret(Default::default())?;
//! assert_eq!(interpretation.structure().atom_count(), 3);
//!
//! let (written, report) = mmcif::write([interpretation.structure()], Default::default())?;
//! assert!(report.omissions().is_empty());
//! let reread = mmcif::parse_str(&written, Default::default())?.interpret(Default::default())?;
//! assert_eq!(reread.structure().atom_count(), 3);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub use crate::io::mmcif_entry_points::{parse_reader, parse_str, write, write_to};
pub use crate::io::{
    MmcifAtomField, MmcifBlock, MmcifBlockSource, MmcifCellReason, MmcifDocument,
    MmcifEntityClassifications, MmcifEntry, MmcifGeometryField, MmcifItem, MmcifLoopTable,
    MmcifParseOptions, MmcifValue, MmcifWriteIssue, MmcifWriteOptions,
};
