//! The `fixed_u32_id!` macro that defines the crate's `u32` identifier types.

macro_rules! fixed_u32_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(u32);

        impl $name {
            /// Wraps a raw value. Nothing is validated here; a lookup with a
            /// value its owner does not hold returns `None`.
            pub const fn new(raw: u32) -> Self {
                Self(raw)
            }

            /// Returns the raw value. It is this identifier's stable
            /// encoding for persistence: `new(id.raw())` restores `id`, and
            /// the encoding does not change between Kekule versions.
            pub const fn raw(self) -> u32 {
                self.0
            }

            /// Returns the raw value as a `usize` position.
            pub const fn index(self) -> usize {
                self.0 as usize
            }
        }
    };
    ($(#[$meta:meta])* $name:ident, $display:literal) => {
        fixed_u32_id!($(#[$meta])* $name);

        impl std::fmt::Display for $name {
            fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(formatter, concat!($display, " {}"), self.0)
            }
        }
    };
}
