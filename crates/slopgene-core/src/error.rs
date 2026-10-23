//! The error type shared by every SlopGene operation.

/// Everything that can go wrong while reading, editing or analysing a sequence.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// A file or string could not be understood.
    #[error("parse error: {0}")]
    Parse(String),

    /// An argument was out of range or otherwise unusable (bad coordinates, bad bases, ...).
    #[error("invalid argument: {0}")]
    Invalid(String),

    /// A named thing (enzyme, feature, primer, file format) does not exist.
    #[error("not found: {0}")]
    NotFound(String),

    /// The input is valid but uses something SlopGene does not handle yet.
    #[error("unsupported: {0}")]
    Unsupported(String),
}

pub type Result<T> = std::result::Result<T, Error>;

pub(crate) fn parse_err(msg: impl Into<String>) -> Error {
    Error::Parse(msg.into())
}

pub(crate) fn invalid(msg: impl Into<String>) -> Error {
    Error::Invalid(msg.into())
}
