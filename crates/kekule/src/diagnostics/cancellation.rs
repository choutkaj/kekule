use std::fmt;
use std::sync::Arc;

/// A host's request to stop a cancellable operation partway.
///
/// It wraps a check that the operation polls between its stages and
/// periodically inside its long loops, a few times per millisecond of work at
/// most. Once the check returns `true`, the operation stops at its next poll
/// and fails with [`Cancelled`], returning no partial result. Operations
/// without a cancellation never poll. The check runs on the thread doing the
/// work, so it should be cheap, such as reading an atomic flag.
///
/// Clones share one check. Equality is identity: clones are equal, while two
/// separately constructed cancellations are not, whatever their checks do.
///
/// ```
/// use std::sync::atomic::{AtomicBool, Ordering};
/// use std::sync::Arc;
/// use kekule::diagnostics::Cancellation;
///
/// let flag = Arc::new(AtomicBool::new(false));
/// let cancellation = Cancellation::from_flag(Arc::clone(&flag));
/// assert!(!cancellation.is_cancelled());
/// flag.store(true, Ordering::Relaxed);
/// assert!(cancellation.is_cancelled());
/// ```
#[derive(Clone)]
pub struct Cancellation {
    check: Arc<dyn Fn() -> bool + Send + Sync>,
}

impl Cancellation {
    /// A cancellation requested once `check` returns `true`.
    pub fn new(check: impl Fn() -> bool + Send + Sync + 'static) -> Self {
        Self {
            check: Arc::new(check),
        }
    }

    /// A cancellation requested once `flag` is set.
    pub fn from_flag(flag: Arc<std::sync::atomic::AtomicBool>) -> Self {
        Self::new(move || flag.load(std::sync::atomic::Ordering::Relaxed))
    }

    /// Whether cancellation has been requested, running the check once.
    pub fn is_cancelled(&self) -> bool {
        (self.check)()
    }

    /// Fails with [`Cancelled`] once cancellation has been requested, running
    /// the check once; for cancellable operations outside this crate.
    ///
    /// # Errors
    ///
    /// Returns [`Cancelled`] when the check returns `true`.
    pub fn check(&self) -> Result<(), Cancelled> {
        if self.is_cancelled() {
            Err(Cancelled)
        } else {
            Ok(())
        }
    }

    /// Polls `cancellation`, which is `None` for an operation that cannot be
    /// cancelled.
    #[inline]
    pub(crate) fn poll(cancellation: Option<&Self>) -> Result<(), Cancelled> {
        match cancellation {
            Some(cancellation) if cancellation.is_cancelled() => Err(Cancelled),
            _ => Ok(()),
        }
    }

    /// Polls `cancellation` on every [`Self::POLL_STRIDE`]th step of a loop
    /// whose steps are too cheap to poll each.
    #[inline]
    pub(crate) fn poll_every(cancellation: Option<&Self>, step: usize) -> Result<(), Cancelled> {
        if cancellation.is_some() && step.is_multiple_of(Self::POLL_STRIDE) {
            Self::poll(cancellation)
        } else {
            Ok(())
        }
    }

    /// Steps between polls of a loop with cheap steps.
    pub(crate) const POLL_STRIDE: usize = 1024;
}

impl fmt::Debug for Cancellation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("Cancellation")
    }
}

impl PartialEq for Cancellation {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.check, &other.check)
    }
}

impl Eq for Cancellation {}

// Kekule only calls the check and reads its answer, so a panic inside it
// cannot leave Kekule state broken, and options holding a cancellation stay
// unwind safe as they were without one.
impl std::panic::UnwindSafe for Cancellation {}
impl std::panic::RefUnwindSafe for Cancellation {}

/// An operation stopped because its [`Cancellation`] was requested.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub struct Cancelled;

impl fmt::Display for Cancelled {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("The operation was cancelled.")
    }
}

impl std::error::Error for Cancelled {}
