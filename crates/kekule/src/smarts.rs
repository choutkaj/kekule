//! SMARTS parsing into syntax-independent substructure queries.
//!
//! [`smarts::parse_str`](crate::smarts::parse_str) returns a [`query::QueryGraph`](crate::query::QueryGraph), the canonical query
//! representation that [`substructure`](crate::substructure) matching consumes. Parsing a query and
//! matching it are separate operations; matching does not perceive the target
//! implicitly. The dialect, including its explicit exclusions, is documented on
//! [`query`](crate::query).
//!
//! ```
//! use kekule::{smarts, smiles, substructure::PreparedTarget};
//!
//! let mut succinic_acid = smiles::to_molecules("OC(=O)CCC(=O)O")?.remove(0);
//! succinic_acid.perceive()?;
//! let carboxylic_acid = smarts::parse_str("C(=O)[OX2H1]", Default::default())?;
//! let target = PreparedTarget::new(&succinic_acid);
//! let matches = target.all(&carboxylic_acid, Default::default())?;
//! assert_eq!(matches.len(), 2);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub use crate::query::smarts::{parse_str, SmartsParseOptions};
