//! Fast, native reader for Juicebox `.hic` files.
//!
//! The reader performs indexed random access. Local files use `pread` and blocks
//! are decompressed in parallel; HTTP(S) inputs use byte range requests.

mod block;
mod error;
mod format;
mod io;
mod v10;

pub use error::{Error, Result};
pub use format::{
    Batch, Chromosome, ContactRecord, HicFile, MatrixType, Normalization, NormalizationEntry,
    PreparedQuery, RawContactRecord, RawValue, Unit,
};
/// Read a sparse matrix window into memory.
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

/// Stream a sparse matrix window to a callback in deterministic block order.
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

/// Read a matrix window into a row-major dense matrix.
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

/// Read a stored normalization vector. `norm` must not be `NONE`.
pub fn normalization_vector(
    file: &str,
    chromosome: &str,
    unit: Unit,
    resolution: i32,
    norm: &Normalization,
) -> Result<Vec<f64>> {
    HicFile::open(file)?.normalization_vector(chromosome, unit, resolution, norm)
}

/// Read a stored expected-value vector. `norm` may be `NONE`.
pub fn expected_vector(
    file: &str,
    chromosome: &str,
    unit: Unit,
    resolution: i32,
    norm: &Normalization,
) -> Result<Vec<f64>> {
    HicFile::open(file)?.expected_vector(chromosome, unit, resolution, norm)
}

/// Read exact V10 raw records, preserving integer counts and float scores
/// separately. Only supported for V10 files.
pub fn raw_records(
    file: &str,
    chr1: &str,
    chr2: &str,
    unit: Unit,
    resolution: i32,
) -> Result<Vec<RawContactRecord>> {
    HicFile::open(file)?.raw_records(chr1, chr2, unit, resolution)
}
