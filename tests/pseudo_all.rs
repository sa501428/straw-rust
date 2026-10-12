use straw::{HicFile, RawValue, Unit};
#[test]
fn scaled_all_layout_and_exact_contacts() {
    let path = format!("{}/tests/data/pseudo-all.hic", env!("CARGO_MANIFEST_DIR"));
    let hic = HicFile::open(path).unwrap();
    let all = hic.all_overview().unwrap().unwrap();
    assert_eq!(all.chromosome_index, 2);
    assert_eq!(all.finest_resolution, 1);
    assert_eq!(all.effective_resolution, 10);
    assert_eq!(all.bin_count, 18);
    assert_eq!(all.source_offsets, vec![0, 10]);
    assert_eq!(all.chromosome_bins, vec![10, 8]);
    let records = hic.raw_records("All", "All", Unit::BP, 1).unwrap();
    assert_eq!(records.len(), 3);
    assert_eq!(records[0].value, RawValue::Count((1u64 << 53) + 17));
    let collapsed = hic.raw_records("All", "All", Unit::BP, 100).unwrap();
    assert_eq!(collapsed.len(), 1);
    assert_eq!(collapsed[0].value, RawValue::Count((1u64 << 53) + 49));
    assert!(hic
        .raw_records("All", "chrA", Unit::BP, 1)
        .unwrap()
        .is_empty());
}
