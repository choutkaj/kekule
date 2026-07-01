mod hierarchy;

// Compatibility/convenience re-export; prefer `molecules::mmcif::read_str` for new code.
pub use crate::io::{read_mmcif_str, MmcifParseError, MmcifParseOptions};
pub use hierarchy::*;
