//! Samples grouped by a host's grouping, such as the bins of a histogram.

/// The sample positions of each group, stored flat: every group's members,
/// group after group, by increasing position.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Grouping {
    /// Where each group's members start in `members`, and one past the last.
    starts: Vec<usize>,
    members: Vec<usize>,
}

impl Grouping {
    /// The groups of `groups` that `group_of` puts each sample position in,
    /// leaving out the samples in no group.
    ///
    /// # Errors
    ///
    /// The first group of `group_of` that is not below `groups`.
    pub(crate) fn new(group_of: &[Option<u32>], groups: usize) -> Result<Self, usize> {
        let mut starts = vec![0_usize; groups + 1];
        for group in group_of.iter().flatten() {
            let group = *group as usize;
            if group >= groups {
                return Err(group);
            }
            starts[group + 1] += 1;
        }
        for group in 0..groups {
            starts[group + 1] += starts[group];
        }
        let mut members = vec![0_usize; starts[groups]];
        let mut next = starts[..groups].to_vec();
        for (position, group) in group_of.iter().enumerate() {
            if let Some(group) = group {
                let slot = &mut next[*group as usize];
                members[*slot] = position;
                *slot += 1;
            }
        }
        Ok(Self { starts, members })
    }

    /// The number of groups.
    pub(crate) fn len(&self) -> usize {
        self.starts.len() - 1
    }

    /// The sample positions of group `group`, by increasing position.
    pub(crate) fn members(&self, group: usize) -> &[usize] {
        &self.members[self.starts[group]..self.starts[group + 1]]
    }
}
