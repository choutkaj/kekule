//! XYZ parsing and structure interpretation.
//!
//! An XYZ file is a sequence of frames, each an atom count, a comment line,
//! and one element and three coordinates per atom. It states no bonds, so
//! interpretation infers connectivity from interatomic distances with the
//! rule macromolecular interpretation applies to residues without a
//! template: two atoms are bonded when they lie within the sum of their
//! covalent radii (Cordero et al. 2008) plus 0.45 Å, and farther apart than
//! 0.4 Å; an element without a tabulated radius is taken as 0.77 Å.
//!
//! Interpretation never guesses chemistry the file does not state. Every
//! inferred bond is single, except that a bond between a metal and a
//! non-metal is dative, as for PDB `CONECT` records; bond orders, formal
//! charges, radicals, and stereo are not assigned. XYZ lists every atom, so
//! every atom has a fixed implicit hydrogen count of zero and perception
//! never adds hydrogens. Each connected component becomes one molecule
//! instance, under a synthetic hierarchy of one `UNL` residue per molecule
//! in chain `A`. The [`XyzInterpretationReport`]
//! lists every inferred bond with its distance and cutoff, atoms that used
//! the fallback radius, atoms too close to bond, and omitted content. The
//! comment line is kept verbatim and not interpreted, and columns after the
//! coordinates, such as extended XYZ properties, are rejected unless
//! [`XyzInterpretOptions::with_extra_columns`]
//! reports them as omitted.
//!
//! One frame becomes a [`Structure`](crate::structure::Structure); several frames of
//! the same atoms become a [`StructureSet`](crate::structure::StructureSet) without
//! an order, and every frame must infer the same bonds. Interpretation keeps atom lines in source order as dense atom
//! order, so a trajectory reader such as [`crate::traj::io::xyz::XyzReader`] can
//! stream the frames of the same file against the interpreted topology when
//! the file spells every element as a symbol in its standard letter case,
//! the only element tokens that reader accepts; interpretation also accepts
//! other letter cases and atomic numbers.
//!
//! ```
//! use kekule::molecule::BondOrder;
//! use kekule::xyz;
//!
//! let text = "3\nwater\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n";
//! let document = xyz::parse_str(text, Default::default())?;
//! let interpretation = document.interpret(Default::default())?;
//! assert_eq!(interpretation.topology().instance_count(), 1);
//! let report = interpretation.report();
//! assert_eq!(report.comment(), "water");
//! assert_eq!(report.inferred_bonds().len(), 2);
//! assert!(report
//!     .inferred_bonds()
//!     .iter()
//!     .all(|bond| bond.order() == BondOrder::Single));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub use crate::io::xyz_entry_points::{parse_reader, parse_str};
pub use crate::io::{
    XyzAtom, XyzDocument, XyzFrame, XyzFrameSelection, XyzInferredBond, XyzInterpretIssue,
    XyzInterpretOptions, XyzInterpretation, XyzInterpretationReport, XyzParseOptions,
    XyzStructureSetInterpretOptions, XyzStructureSetInterpretation, DEFAULT_XYZ_MAX_ATOM_CONTACTS,
};
