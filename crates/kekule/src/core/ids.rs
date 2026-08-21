#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct FixedIdCapacityError;

pub(crate) fn checked_raw_id(index: usize) -> Result<u32, FixedIdCapacityError> {
    u32::try_from(index).map_err(|_| FixedIdCapacityError)
}

pub(crate) fn checked_fixed_id_collection_len(
    current: usize,
    additional: usize,
) -> Result<(), FixedIdCapacityError> {
    let final_len = current
        .checked_add(additional)
        .ok_or(FixedIdCapacityError)?;
    if final_len == 0 {
        return Ok(());
    }
    checked_raw_id(final_len - 1).map(|_| ())
}

fixed_u32_id!(AtomId, "a");
fixed_u32_id!(BondId, "b");
fixed_u32_id!(StereoElementId, "s");
fixed_u32_id!(StereoGroupId, "sg");
