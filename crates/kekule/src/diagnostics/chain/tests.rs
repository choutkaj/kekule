use super::*;
use crate::diagnostics::{Format, FormatError};

#[test]
fn error_chains_read_as_sentences() {
    #[derive(Debug)]
    struct Outer(std::io::Error);
    impl fmt::Display for Outer {
        fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            formatter.write_str("Cannot build the solvated structure.")
        }
    }
    impl std::error::Error for Outer {
        fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
            Some(&self.0)
        }
    }

    let error = Outer(std::io::Error::other("disk full"));
    assert_eq!(
        ErrorChain::new(&error).to_string(),
        "Cannot build the solvated structure. Disk full."
    );
    assert_eq!(
        format!("{:#}", ErrorChain::new(&error)),
        "Cannot build the solvated structure.
Disk full."
    );

    let sink = FormatError::sink(
        Format::Pdb,
        std::io::Error::other("Access is denied. (os error 5)"),
    );
    assert_eq!(
        ErrorChain::new(&sink).to_string(),
        "Cannot write the PDB output. Access is denied. (os error 5)."
    );
}
