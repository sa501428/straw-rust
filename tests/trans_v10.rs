use straw::{HicFile, MatrixType, Normalization, Unit};

// Independent byte-level fixtures from straw/C++/tests/test_v10.py:
// fixture(path, trans=True, frag=True), plus score=True and f32 bits for (1, 1, 5).
// chrA=80 bp, chrB=70 bp; sum=7, source=10 bp, derived=20 bp.
// No cis expected indexes. Normalization divisors are A=2/4, B=3/6.
#[test]
fn trans_oe_uses_whole_matrix_mean_without_expected_vectors() {
    for name in ["trans-v10.hic", "trans-score-v10.hic"] {
        let path = format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"));
        let hic = HicFile::open(path).unwrap();
        let frag = hic
            .records(
                MatrixType::Oe,
                Normalization::none(),
                "chrB",
                "chrA",
                Unit::Frag,
                1,
            )
            .unwrap();
        assert_eq!(frag.len(), 3);
        assert!((f64::from(frag[0].counts) - 64.0 / 7.0).abs() < 1e-5);
        for (resolution, mean, product) in [(10, 7.0 / 56.0, 6.0), (20, 7.0 / 16.0, 24.0)] {
            for norm in [Normalization::none(), Normalization::new("VC")] {
                for mt in [MatrixType::Oe, MatrixType::Expected] {
                    for (a, b) in [("chrA", "chrB"), ("chrB", "chrA")] {
                        let records = hic
                            .records(mt, norm.clone(), a, b, Unit::BP, resolution)
                            .unwrap();
                        assert_eq!(records.len(), 3);
                        for (record, count) in records.iter().zip([1.0, 1.0, 5.0]) {
                            let expected = if mt == MatrixType::Expected {
                                mean
                            } else {
                                count / mean / if norm.is_none() { 1.0 } else { product }
                            };
                            assert!((f64::from(record.counts) - expected).abs() < 1e-5);
                        }
                        let forward = if resolution == 10 {
                            vec![(0, 0), (10, 20), (30, 30)]
                        } else {
                            vec![(0, 0), (0, 20), (20, 20)]
                        };
                        let coordinates: Vec<_> =
                            records.iter().map(|r| (r.bin_x, r.bin_y)).collect();
                        assert_eq!(
                            coordinates,
                            if a == "chrA" {
                                forward
                            } else {
                                forward.into_iter().map(|(x, y)| (y, x)).collect()
                            }
                        );
                    }
                }
            }
            let prepared = hic
                .prepare(
                    MatrixType::Oe,
                    Normalization::none(),
                    "chrB",
                    "chrA",
                    Unit::BP,
                    resolution,
                )
                .unwrap();
            let window = prepared.window([0, 10, 0, 10]).unwrap();
            assert_eq!(window.len(), 1);
            assert!((f64::from(window[0].counts) - 1.0 / mean).abs() < 1e-5);
            let matrix = hic
                .matrix(
                    MatrixType::Oe,
                    Normalization::none(),
                    "chrB:0:40",
                    "chrA:0:40",
                    Unit::BP,
                    resolution,
                )
                .unwrap();
            assert!((f64::from(matrix[0][0]) - 1.0 / mean).abs() < 1e-5);
        }
    }
}

#[test]
fn zero_and_nonfinite_trans_means_are_unavailable() {
    let original = include_bytes!("data/trans-score-v10.hic");
    let footer = u64::from_le_bytes(original[16..24].try_into().unwrap()) as usize;
    let matrix =
        u64::from_le_bytes(original[footer + 32..footer + 40].try_into().unwrap()) as usize;
    let path = std::env::temp_dir().join(format!("straw-trans-mean-{}.hic", std::process::id()));
    for mean in [0.0f64, f64::NAN, f64::INFINITY] {
        let mut data = original.to_vec();
        data[matrix + 44..matrix + 52].copy_from_slice(&mean.to_bits().to_le_bytes());
        std::fs::write(&path, data).unwrap();
        let hic = HicFile::open(path.to_str().unwrap()).unwrap();
        for mt in [MatrixType::Oe, MatrixType::Expected] {
            assert!(hic
                .records(mt, Normalization::none(), "chrA", "chrB", Unit::BP, 10)
                .unwrap()
                .is_empty());
        }
    }
    std::fs::remove_file(path).unwrap();
}
