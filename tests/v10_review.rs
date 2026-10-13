use straw::{HicFile, MatrixType, Normalization, RawValue, Unit};
fn file(name: &str) -> HicFile {
    HicFile::open(format!("{}/tests/data/{name}.hic", env!("CARGO_MANIFEST_DIR"))).unwrap()
}
#[test]
fn inclusive_sparse_point_and_dense_bounds() {
    let hic = file("trans-v10");
    for resolution in [10, 20] {
        let x = if resolution == 10 { 10 } else { 0 };
        let a = format!("chrA:{x}:{x}");
        let b = "chrB:20:20";
        let records = hic.records(MatrixType::Observed, Normalization::none(), &a, b, Unit::BP, resolution).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!((records[0].bin_x, records[0].bin_y, records[0].counts), (x, 20, 1.));
        let raw = hic.raw_records(&a, b, Unit::BP, resolution).unwrap();
        assert_eq!(raw.len(), 1);
        let matrix = hic.matrix(MatrixType::Observed, Normalization::none(), &a, b, Unit::BP, resolution).unwrap();
        assert_eq!(matrix, vec![vec![1.]]);
    }
    assert!(hic.records(MatrixType::Observed, Normalization::none(), "chrA:11:19", "chrB:20:20", Unit::BP, 10).unwrap().is_empty());
    let matrix = hic.matrix(MatrixType::Observed, Normalization::none(), "chrA:1:10", "chrB:11:20", Unit::BP, 10).unwrap();
    assert_eq!(matrix, vec![vec![0., 0.], vec![0., 1.]]);
    let frag = hic.matrix(MatrixType::Observed, Normalization::none(), "chrA", "chrB", Unit::Frag, 1).unwrap();
    assert_eq!(frag.len(), 9); // eight fragments; inclusive numeric end produces one empty bin
}
#[test]
fn derived_buffers_subdivide_without_losing_or_duplicating_cells() {
    for name in ["large-derived-count", "large-derived-score"] {
        let hic = file(name);
        let dense = hic.matrix(MatrixType::Observed, Normalization::none(), "chrA:1:2", "chrA:0:2", Unit::BP, 2).unwrap();
        assert_eq!(dense, vec![vec![0., 0.], vec![1., 0.]]);
        let records = hic.raw_records("chrA", "chrA", Unit::BP, 2).unwrap();
        assert_eq!(records.len(), 70000);
        for (i, r) in records.iter().enumerate() {
            assert_eq!((r.bin_x, r.bin_y), (i as u64, i as u64 + 1));
            match r.value { RawValue::Count(v) => assert_eq!(v, 1), RawValue::Score(v) => assert_eq!(v, 1.) }
        }
        let overlap = hic.records(MatrixType::Observed, Normalization::none(), "chrA:2:139998", "chrA:0:139996", Unit::BP, 2).unwrap();
        assert_eq!(overlap.len(), 69999);
        let actual: std::collections::BTreeSet<_> = overlap.iter().map(|r| (r.bin_x, r.bin_y)).collect();
        let mut expected: std::collections::BTreeSet<_> = (1..69998).map(|i| (2*i, 2*i+2)).collect();
        expected.insert((2,0)); expected.insert((139998,139996));
        assert_eq!(actual, expected);
        assert_eq!(overlap[0].bin_x, 2);
        assert_eq!(overlap.last().unwrap().bin_x, 139998);
    }
    let hic = file("large-score-fanin");
    let records = hic.raw_records("chrA", "chrB", Unit::BP, 1000).unwrap();
    assert_eq!(records.len(), 1);
    match records[0].value { RawValue::Score(v) => assert_eq!(v, 70000.), _ => panic!("score expected") }
}
