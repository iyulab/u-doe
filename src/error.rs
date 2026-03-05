/// Errors produced by u-doe operations.
#[derive(Debug, Clone, PartialEq)]
pub enum DoeError {
    /// Invalid number of factors for this design type.
    InvalidFactorCount { min: usize, max: usize, got: usize },
    /// Insufficient response data.
    InsufficientResponses { expected: usize, got: usize },
    /// Specification limits are invalid.
    InvalidSpecification(String),
    /// Matrix operation failed.
    MatrixError(String),
    /// Requested design not supported in v0.1.
    UnsupportedDesign(String),
}

impl std::fmt::Display for DoeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DoeError::InvalidFactorCount { min, max, got } => write!(
                f, "invalid factor count: expected {min}..={max}, got {got}"
            ),
            DoeError::InsufficientResponses { expected, got } => write!(
                f, "insufficient responses: expected {expected}, got {got}"
            ),
            DoeError::InvalidSpecification(msg) => write!(f, "invalid specification: {msg}"),
            DoeError::MatrixError(msg) => write!(f, "matrix error: {msg}"),
            DoeError::UnsupportedDesign(msg) => write!(f, "unsupported design: {msg}"),
        }
    }
}

impl std::error::Error for DoeError {}
