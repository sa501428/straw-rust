//! Regression tests for the user-requested-axis contact orientation contract:
//! the first returned coordinate always corresponds to the first chromosome
//! or region the caller passed, regardless of on-disk chromosome-index order.

use std::collections::HashMap;
use straw::{straw, ContactRecord, MatrixType, Normalization, Unit};

fn fixture() -> String {
    concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../straw/R/inst/extdata/test.hic"
    )
    .to_string()
}

fn multiset(records: Vec<ContactRecord>) -> HashMap<(i32, i32, u32), usize> {
    let mut map = HashMap::new();
    for r in records {
        *map.entry((r.bin_x, r.bin_y, r.counts.to_bits()))
            .or_insert(0) += 1;
    }
    map
}

fn transposed(map: &HashMap<(i32, i32, u32), usize>) -> HashMap<(i32, i32, u32), usize> {
    map.iter()
        .map(|(&(x, y, v), &count)| ((y, x, v), count))
        .collect()
}

#[test]
fn interchromosomal_reversal_transposes() {
    let path = fixture();
    let forward = straw(
        MatrixType::Observed,
        Normalization::none(),
        &path,
        "1",
        "2",
        Unit::BP,
        2_500_000,
    )
    .unwrap();
    let reverse = straw(
        MatrixType::Observed,
        Normalization::none(),
        &path,
        "2",
        "1",
        Unit::BP,
        2_500_000,
    )
    .unwrap();
    assert!(
        !forward.is_empty(),
        "fixture returned no interchromosomal contacts"
    );
    assert_eq!(multiset(reverse), transposed(&multiset(forward)));
}

#[test]
fn cis_asymmetric_window_reflects() {
    let path = fixture();
    let upper = straw(
        MatrixType::Observed,
        Normalization::none(),
        &path,
        "1:0:5000000",
        "1:10000000:15000000",
        Unit::BP,
        2_500_000,
    )
    .unwrap();
    let lower = straw(
        MatrixType::Observed,
        Normalization::none(),
        &path,
        "1:10000000:15000000",
        "1:0:5000000",
        Unit::BP,
        2_500_000,
    )
    .unwrap();
    assert!(
        !upper.is_empty(),
        "fixture returned no asymmetric cis contacts"
    );
    assert_eq!(multiset(lower), transposed(&multiset(upper)));
}

#[test]
fn dense_matrix_reversal_is_transpose() {
    let path = fixture();
    let forward = straw::straw_as_matrix(
        MatrixType::Observed,
        Normalization::none(),
        &path,
        "1",
        "2",
        Unit::BP,
        2_500_000,
    )
    .unwrap();
    let reverse = straw::straw_as_matrix(
        MatrixType::Observed,
        Normalization::none(),
        &path,
        "2",
        "1",
        Unit::BP,
        2_500_000,
    )
    .unwrap();
    assert_eq!(forward.len(), reverse[0].len());
    assert_eq!(forward[0].len(), reverse.len());
    for (r, row) in forward.iter().enumerate() {
        for (c, &value) in row.iter().enumerate() {
            assert_eq!(value, reverse[c][r], "mismatch at ({r}, {c})");
        }
    }
}
