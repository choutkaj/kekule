fixed_u32_id!(
    /// A molecule-local atom ID, owned by one [`Molecule`](super::Molecule)
    /// or editor draft.
    ///
    /// In a published molecule it is also the dense atom index `0..atom_count`;
    /// draft IDs stay stable but may be sparse until publication renumbers them.
    AtomId,
    "atom"
);
fixed_u32_id!(
    /// A molecule-local bond ID, owned by one [`Molecule`](super::Molecule)
    /// or editor draft.
    ///
    /// In a published molecule it is also the dense bond index `0..bond_count`;
    /// draft IDs stay stable but may be sparse until publication renumbers them.
    BondId,
    "bond"
);
fixed_u32_id!(
    /// A molecule-local stereo element ID, owned by one
    /// [`Molecule`](super::Molecule) or editor draft.
    ///
    /// Published IDs are dense; draft IDs stay stable but may be sparse until
    /// publication renumbers them.
    StereoElementId,
    "stereo element"
);
fixed_u32_id!(
    /// A molecule-local enhanced stereo group ID, owned by one
    /// [`Molecule`](super::Molecule) or editor draft.
    ///
    /// Published IDs are dense; draft IDs stay stable but may be sparse until
    /// publication renumbers them.
    StereoGroupId,
    "stereo group"
);
