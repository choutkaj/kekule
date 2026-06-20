use crate::*;

pub(crate) fn boxed_error(message: impl Into<String>) -> Box<dyn Error> {
    std::io::Error::other(message.into()).into()
}
