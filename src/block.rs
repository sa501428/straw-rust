use crate::io::{Bytes, RandomAccess};
use crate::{Error, Result};
use flate2::read::ZlibDecoder;
use std::io::{Cursor, Read};
use std::sync::Arc;

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct IndexEntry {
    pub position: u64,
    pub size: u64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RawRecord {
    pub bin_x: i32,
    pub bin_y: i32,
    pub counts: f32,
}

pub(crate) fn read_block(
    source: &Arc<dyn RandomAccess>,
    entry: IndexEntry,
    version: i32,
) -> Result<Vec<RawRecord>> {
    let decoded = decompress(source, entry)?;
    decode_block(&decoded, version)
}

fn decompress(source: &Arc<dyn RandomAccess>, entry: IndexEntry) -> Result<Vec<u8>> {
    if entry.size == 0 {
        return Ok(Vec::new());
    }
    let compressed = source.read_exact_at(
        entry.position,
        usize::try_from(entry.size).map_err(|_| Error::CorruptBlock("block too large".into()))?,
    )?;
    let mut decoded = Vec::with_capacity(compressed.len().saturating_mul(3));
    if compressed.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]) {
        let mut decoder = ruzstd::StreamingDecoder::new(Cursor::new(compressed))
            .map_err(|e| Error::CorruptBlock(format!("zstd: {e:?}")))?;
        decoder.read_to_end(&mut decoded)?;
    } else {
        ZlibDecoder::new(compressed.as_slice()).read_to_end(&mut decoded)?;
    }
    Ok(decoded)
}

pub(crate) fn record_count(
    source: &Arc<dyn RandomAccess>,
    entry: IndexEntry,
    version: i32,
) -> Result<u64> {
    let decoded = decompress(source, entry)?;
    if decoded.len() < 4 {
        return Err(Error::CorruptBlock("missing record count".into()));
    }
    let count = i32::from_le_bytes(decoded[..4].try_into().unwrap());
    if count < 0 {
        return Err(Error::CorruptBlock("negative record count".into()));
    }
    let _ = version;
    Ok(count as u64)
}

fn decode_block(data: &[u8], version: i32) -> Result<Vec<RawRecord>> {
    let mut r = Bytes::new(data);
    let declared = r.i32()?;
    if declared < 0 {
        return Err(Error::CorruptBlock("negative record count".into()));
    }
    let mut out = Vec::with_capacity(declared as usize);
    if version < 7 {
        for _ in 0..declared {
            out.push(RawRecord {
                bin_x: r.i32()?,
                bin_y: r.i32()?,
                counts: r.f32()?,
            });
        }
        return Ok(out);
    }
    let x_offset = r.i32()?;
    let y_offset = r.i32()?;
    let short_counts = r.u8()? == 0;
    let (short_x, short_y) = if version > 8 {
        (r.u8()? == 0, r.u8()? == 0)
    } else {
        (true, true)
    };
    let encoding = r.u8()?;
    let (all_one, delta_columns) = if version > 9 {
        (r.u8()? != 0, r.u8()? != 0)
    } else {
        (false, false)
    };

    match encoding {
        1 => {
            let rows = if short_y { r.i16()? as i32 } else { r.i32()? };
            if rows < 0 {
                return Err(Error::CorruptBlock("negative row count".into()));
            }
            for _ in 0..rows {
                let dy = if short_y { r.i16()? as i32 } else { r.i32()? };
                let cols = if short_x { r.i16()? as i32 } else { r.i32()? };
                if cols < 0 {
                    return Err(Error::CorruptBlock("negative column count".into()));
                }
                let mut previous = 0i32;
                for _ in 0..cols {
                    let dx = if short_x { r.i16()? as i32 } else { r.i32()? };
                    previous = if delta_columns {
                        previous
                            .checked_add(dx)
                            .ok_or_else(|| Error::CorruptBlock("column overflow".into()))?
                    } else {
                        dx
                    };
                    out.push(RawRecord {
                        bin_x: x_offset + previous,
                        bin_y: y_offset + dy,
                        counts: read_count(&mut r, all_one, short_counts)?,
                    });
                }
            }
        }
        2 => {
            let points = r.i32()?;
            let width = r.i16()? as i32;
            if points < 0 || width <= 0 {
                return Err(Error::CorruptBlock("invalid dense dimensions".into()));
            }
            for i in 0..points {
                let x = x_offset + i % width;
                let y = y_offset + i / width;
                if short_counts {
                    let c = r.i16()?;
                    if c != i16::MIN {
                        out.push(RawRecord {
                            bin_x: x,
                            bin_y: y,
                            counts: c as f32,
                        });
                    }
                } else {
                    let c = r.f32()?;
                    if !c.is_nan() {
                        out.push(RawRecord {
                            bin_x: x,
                            bin_y: y,
                            counts: c,
                        });
                    }
                }
            }
        }
        other => return Err(Error::CorruptBlock(format!("unknown encoding {other}"))),
    }
    // Some dense blocks declare grid slots, not emitted records.
    out.shrink_to_fit();
    Ok(out)
}

fn read_count(r: &mut Bytes<'_>, all_one: bool, short: bool) -> Result<f32> {
    if all_one {
        Ok(1.0)
    } else if short {
        Ok(r.i16()? as f32)
    } else {
        r.f32()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn version_six_records() {
        let mut b = 1i32.to_le_bytes().to_vec();
        b.extend(2i32.to_le_bytes());
        b.extend(3i32.to_le_bytes());
        b.extend(4.5f32.to_le_bytes());
        assert_eq!(
            decode_block(&b, 6).unwrap(),
            vec![RawRecord {
                bin_x: 2,
                bin_y: 3,
                counts: 4.5
            }]
        );
    }

    #[test]
    fn version_nine_dense_records_skip_sentinels() {
        let mut b = 2i32.to_le_bytes().to_vec();
        b.extend(10i32.to_le_bytes());
        b.extend(20i32.to_le_bytes());
        b.extend([0, 0, 0, 2]);
        b.extend(3i32.to_le_bytes());
        b.extend(3i16.to_le_bytes());
        for count in [5i16, i16::MIN, 7] {
            b.extend(count.to_le_bytes());
        }
        assert_eq!(
            decode_block(&b, 9).unwrap(),
            vec![
                RawRecord {
                    bin_x: 10,
                    bin_y: 20,
                    counts: 5.0
                },
                RawRecord {
                    bin_x: 12,
                    bin_y: 20,
                    counts: 7.0
                },
            ]
        );
    }

    #[test]
    fn version_ten_delta_columns_and_implicit_counts() {
        let mut b = 2i32.to_le_bytes().to_vec();
        b.extend(10i32.to_le_bytes());
        b.extend(20i32.to_le_bytes());
        b.extend([0, 0, 0, 1, 1, 1]);
        b.extend(1i16.to_le_bytes());
        b.extend(2i16.to_le_bytes());
        b.extend(2i16.to_le_bytes());
        b.extend(3i16.to_le_bytes());
        b.extend(4i16.to_le_bytes());
        assert_eq!(
            decode_block(&b, 10).unwrap(),
            vec![
                RawRecord {
                    bin_x: 13,
                    bin_y: 22,
                    counts: 1.0
                },
                RawRecord {
                    bin_x: 17,
                    bin_y: 22,
                    counts: 1.0
                },
            ]
        );
    }
}
