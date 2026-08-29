use crate::format::{HicFile, MatrixType, Normalization, Unit};
use crate::{Error, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ContactFilter {
    #[default]
    All,
    IntraShort,
    IntraLong,
    Inter,
    Intra,
}

#[derive(Clone, Debug)]
pub struct DumpOptions {
    pub matrix_type: MatrixType,
    pub normalization: Normalization,
    pub unit: Unit,
    pub resolution: i32,
    pub compressed: bool,
    pub filter: ContactFilter,
}

/// Write the C++ straw `HICSLICE` format. Records deliberately preserve the
/// C++ struct's padding, making files binary-compatible (20 bytes/record).
pub fn dump(input: &str, output: &str, options: &DumpOptions) -> Result<()> {
    if options.resolution <= 0 {
        return Err(Error::Argument("resolution must be positive".into()));
    }
    let hic = HicFile::open(input)?;
    let chromosomes: Vec<_> = hic
        .chromosomes
        .iter()
        .filter(|c| {
            if hic.version() == 10 {
                c.name != "All" && c.name != "ALL"
            } else {
                c.index > 0
            }
        })
        .cloned()
        .collect();
    let mut keys = BTreeMap::new();
    for (key, chromosome) in chromosomes.iter().enumerate() {
        let key = i16::try_from(key)
            .map_err(|_| Error::Invalid("too many chromosomes for slice format".into()))?;
        keys.insert(chromosome.name.clone(), key);
    }
    let file = File::create(output)?;
    let writer: Box<dyn Write> = if options.compressed {
        Box::new(GzEncoder::new(
            BufWriter::with_capacity(1024 * 1024, file),
            Compression::default(),
        ))
    } else {
        Box::new(BufWriter::with_capacity(1024 * 1024, file))
    };
    write_slice(&hic, &chromosomes, keys, options, writer)
}

fn write_slice(
    hic: &HicFile,
    chromosomes: &[crate::Chromosome],
    keys: BTreeMap<String, i16>,
    options: &DumpOptions,
    mut out: Box<dyn Write>,
) -> Result<()> {
    out.write_all(b"HICSLICE")?;
    out.write_all(&options.resolution.to_le_bytes())?;
    out.write_all(&(i32::try_from(keys.len()).unwrap()).to_le_bytes())?;
    for (name, key) in &keys {
        out.write_all(&(name.len() as i32).to_le_bytes())?;
        out.write_all(name.as_bytes())?;
        out.write_all(&key.to_le_bytes())?;
    }
    for (i, c1) in chromosomes.iter().enumerate() {
        for c2 in chromosomes.iter().skip(i) {
            if options.filter == ContactFilter::Inter && c1.index == c2.index {
                continue;
            }
            if matches!(
                options.filter,
                ContactFilter::Intra | ContactFilter::IntraShort | ContactFilter::IntraLong
            ) && c1.index != c2.index
            {
                continue;
            }
            if hic.version() == 10 {
                let records = hic.records(
                    options.matrix_type,
                    options.normalization.clone(),
                    &c1.name,
                    &c2.name,
                    options.unit,
                    options.resolution,
                )?;
                for rec in records {
                    let bin_x = rec.bin_x / options.resolution;
                    let bin_y = rec.bin_y / options.resolution;
                    if rec.counts <= 0.0
                        || !rec.counts.is_finite()
                        || !keep(
                            bin_x,
                            bin_y,
                            c1.index == c2.index,
                            options.resolution,
                            options.filter,
                        )
                    {
                        continue;
                    }
                    out.write_all(&keys[&c1.name].to_le_bytes())?;
                    out.write_all(&[0, 0])?;
                    out.write_all(&bin_x.to_le_bytes())?;
                    out.write_all(&keys[&c2.name].to_le_bytes())?;
                    out.write_all(&[0, 0])?;
                    out.write_all(&bin_y.to_le_bytes())?;
                    out.write_all(&rec.counts.to_le_bytes())?;
                }
                continue;
            }
            let zoom = match hic.load_zoom(
                c1,
                c2,
                options.matrix_type,
                &options.normalization,
                options.unit,
                options.resolution,
            ) {
                Ok(z) => z,
                Err(Error::MatrixNotFound(_) | Error::ResolutionNotFound { .. }) => continue,
                Err(e) => return Err(e),
            };
            for entry in zoom.blocks.values() {
                for rec in zoom.raw_block(*entry)? {
                    if rec.counts <= 0.0
                        || !rec.counts.is_finite()
                        || !keep(
                            rec.bin_x,
                            rec.bin_y,
                            c1.index == c2.index,
                            options.resolution,
                            options.filter,
                        )
                    {
                        continue;
                    }
                    out.write_all(&keys[&c1.name].to_le_bytes())?;
                    out.write_all(&[0, 0])?;
                    out.write_all(&rec.bin_x.to_le_bytes())?;
                    out.write_all(&keys[&c2.name].to_le_bytes())?;
                    out.write_all(&[0, 0])?;
                    out.write_all(&rec.bin_y.to_le_bytes())?;
                    out.write_all(&rec.counts.to_le_bytes())?;
                }
            }
        }
    }
    out.flush()?;
    Ok(())
}

fn keep(x: i32, y: i32, intra: bool, resolution: i32, filter: ContactFilter) -> bool {
    match filter {
        ContactFilter::All => true,
        ContactFilter::Inter => !intra,
        ContactFilter::Intra => intra,
        ContactFilter::IntraShort => {
            intra && (x - y).unsigned_abs() < (5_000_000 / resolution) as u32
        }
        ContactFilter::IntraLong => {
            intra && (x - y).unsigned_abs() > (5_000_000 / resolution) as u32
        }
    }
}
