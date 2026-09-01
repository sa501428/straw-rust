//! Coverage for the normalization/expected vector, prepared-query batching,
//! and exact-raw-record APIs added to close gaps against the cross-language
//! C API plan (straw/C_API_PLAN.md).

use straw::{ContactRecord, Error, HicFile, MatrixType, Normalization, Unit};

fn fixture() -> String {
    concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../straw/R/inst/extdata/test.hic"
    )
    .to_string()
}

#[test]
fn normalization_vector_matches_records() {
    let hic = HicFile::open(fixture()).unwrap();
    let norm = Normalization::new("KR");
    let vector = hic
        .normalization_vector("1", Unit::BP, 2_500_000, &norm)
        .unwrap();
    assert!(!vector.is_empty());

    // Applying the vector by hand should match a normalized query.
    let raw = hic
        .records(
            MatrixType::Observed,
            Normalization::none(),
            "1",
            "1",
            Unit::BP,
            2_500_000,
        )
        .unwrap();
    let normalized = hic
        .records(MatrixType::Observed, norm, "1", "1", Unit::BP, 2_500_000)
        .unwrap();
    assert_eq!(raw.len(), normalized.len());
    for (r, n) in raw.iter().zip(&normalized) {
        let bx = (r.bin_x / 2_500_000) as usize;
        let by = (r.bin_y / 2_500_000) as usize;
        let expected = r.counts as f64 / (vector[bx] * vector[by]);
        assert!((expected as f32 - n.counts).abs() < 1e-3);
    }
}

#[test]
fn normalization_vector_rejects_none() {
    let hic = HicFile::open(fixture()).unwrap();
    let err = hic
        .normalization_vector("1", Unit::BP, 2_500_000, &Normalization::none())
        .unwrap_err();
    assert!(matches!(err, Error::Argument(_)));
}

#[test]
fn expected_vector_is_nonempty_and_matches_oe() {
    let hic = HicFile::open(fixture()).unwrap();
    let norm = Normalization::none();
    let expected = hic
        .expected_vector("1", Unit::BP, 2_500_000, &norm)
        .unwrap();
    assert!(!expected.is_empty());

    let observed = hic
        .records(
            MatrixType::Observed,
            Normalization::none(),
            "1",
            "1",
            Unit::BP,
            2_500_000,
        )
        .unwrap();
    let oe = hic
        .records(MatrixType::Oe, norm, "1", "1", Unit::BP, 2_500_000)
        .unwrap();
    for (o, e) in observed.iter().zip(&oe) {
        let distance = (o.bin_x.abs_diff(o.bin_y) / 2_500_000) as usize;
        let denom = expected[distance.min(expected.len() - 1)];
        assert!(((o.counts as f64 / denom) as f32 - e.counts).abs() < 1e-3);
    }
}

#[test]
fn prepared_query_window_matches_direct_query() {
    let hic = HicFile::open(fixture()).unwrap();
    let direct = hic
        .records(
            MatrixType::Observed,
            Normalization::none(),
            "1:0:5000000",
            "1:10000000:15000000",
            Unit::BP,
            2_500_000,
        )
        .unwrap();

    let prepared = hic
        .prepare(
            MatrixType::Observed,
            Normalization::none(),
            "1",
            "1",
            Unit::BP,
            2_500_000,
        )
        .unwrap();
    let windowed = prepared
        .window([0, 5_000_000, 10_000_000, 15_000_000])
        .unwrap();

    assert_eq!(multiset(direct), multiset(windowed));
}

#[test]
fn prepared_query_regions_concatenates_in_order() {
    let hic = HicFile::open(fixture()).unwrap();
    let prepared = hic
        .prepare(
            MatrixType::Observed,
            Normalization::none(),
            "1",
            "1",
            Unit::BP,
            2_500_000,
        )
        .unwrap();

    let region_a = [0, 5_000_000, 10_000_000, 15_000_000];
    let region_b = [0, 5_000_000, 0, 5_000_000];
    let batch = prepared.regions(&[region_a, region_b]).unwrap();

    assert_eq!(batch.offsets.len(), 3);
    assert_eq!(batch.offsets[0], 0);
    assert_eq!(batch.offsets[2], batch.x.len());

    let expected_a = prepared.window(region_a).unwrap();
    let expected_b = prepared.window(region_b).unwrap();
    assert_eq!(batch.offsets[1] - batch.offsets[0], expected_a.len());
    assert_eq!(batch.offsets[2] - batch.offsets[1], expected_b.len());
}

#[test]
fn raw_records_are_unsupported_on_legacy_files() {
    let hic = HicFile::open(fixture()).unwrap();
    let err = hic.raw_records("1", "2", Unit::BP, 2_500_000).unwrap_err();
    assert!(matches!(err, Error::Unsupported(_)));
}

fn multiset(records: Vec<ContactRecord>) -> std::collections::HashMap<(i32, i32, u32), usize> {
    let mut map = std::collections::HashMap::new();
    for r in records {
        *map.entry((r.bin_x, r.bin_y, r.counts.to_bits()))
            .or_insert(0) += 1;
    }
    map
}
