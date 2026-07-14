//! Fast, native reader for Juicebox `.hic` files.
//!
//! The reader performs indexed random access. Local files use `pread` and blocks
//! are decompressed in parallel; HTTP(S) inputs use byte range requests.

mod block;
mod error;
mod format;
mod io;
mod slice;

pub use error::{Error, Result};
pub use format::{Chromosome, ContactRecord, HicFile, MatrixType, Normalization, Unit};
pub use slice::{dump, ContactFilter, DumpOptions};

/// Read a sparse matrix slice into memory.
pub fn straw(
    matrix_type: MatrixType,
    norm: Normalization,
    file: &str,
    chr1: &str,
    chr2: &str,
    unit: Unit,
    resolution: i32,
) -> Result<Vec<ContactRecord>> {
    HicFile::open(file)?.records(matrix_type, norm, chr1, chr2, unit, resolution)
}

/// Stream a sparse matrix slice to a callback in deterministic block order.
#[allow(clippy::too_many_arguments)] // Mirrors the established straw API.
pub fn straw_stream<F>(
    matrix_type: MatrixType,
    norm: Normalization,
    file: &str,
    chr1: &str,
    chr2: &str,
    unit: Unit,
    resolution: i32,
    callback: F,
) -> Result<()>
where
    F: FnMut(ContactRecord),
{
    HicFile::open(file)?.stream_records(matrix_type, norm, chr1, chr2, unit, resolution, callback)
}

/// Read a matrix slice into a row-major dense matrix.
pub fn straw_as_matrix(
    matrix_type: MatrixType,
    norm: Normalization,
    file: &str,
    chr1: &str,
    chr2: &str,
    unit: Unit,
    resolution: i32,
) -> Result<Vec<Vec<f32>>> {
    HicFile::open(file)?.matrix(matrix_type, norm, chr1, chr2, unit, resolution)
}

/// Count declared contact records across all chromosome pairs at a BP resolution.
pub fn get_num_records_for_file(file: &str, resolution: i32, inter_only: bool) -> Result<u64> {
    HicFile::open(file)?.count_records(resolution, inter_only)
}

/// Count declared intra-chromosomal records independently for every chromosome.
pub fn get_num_records_for_chromosomes(
    file: &str,
    resolution: i32,
) -> Result<Vec<(Chromosome, u64)>> {
    HicFile::open(file)?.chromosome_record_counts(resolution)
}
