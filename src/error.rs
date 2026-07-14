use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("HTTP error: {0}")]
    Http(#[from] reqwest::Error),
    #[error("invalid hic file: {0}")]
    Invalid(String),
    #[error("unsupported hic version {0}; versions 6 and newer are supported")]
    UnsupportedVersion(i32),
    #[error("chromosome {0} not found in the file")]
    ChromosomeNotFound(String),
    #[error("matrix {0} not found")]
    MatrixNotFound(String),
    #[error("resolution {resolution} {unit} not found")]
    ResolutionNotFound { resolution: i32, unit: String },
    #[error(
        "normalization vector {norm} for chromosome {chromosome} at {resolution} {unit} not found"
    )]
    NormalizationNotFound {
        norm: String,
        chromosome: i32,
        resolution: i32,
        unit: String,
    },
    #[error("expected-value vector at {resolution} {unit} not found")]
    ExpectedNotFound { resolution: i32, unit: String },
    #[error("invalid argument: {0}")]
    Argument(String),
    #[error("corrupt block: {0}")]
    CorruptBlock(String),
}

pub type Result<T> = std::result::Result<T, Error>;
