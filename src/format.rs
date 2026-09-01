use crate::block::{read_block, record_count, IndexEntry, RawRecord};
use crate::io::{open_source, Bytes, RandomAccess, Reader};
use crate::v10::V10File;
use crate::{Error, Result};
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt;
use std::str::FromStr;
use std::sync::{Arc, OnceLock};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Chromosome {
    pub name: String,
    pub index: i32,
    pub length: i64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactRecord {
    pub bin_x: i32,
    pub bin_y: i32,
    pub counts: f32,
}

/// A V10 raw value, preserving whether the file stores an exact integer count
/// or a legacy-precision float score for this bin pair. Only V10 files
/// distinguish the two; see [`HicFile::raw_records`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RawValue {
    Count(u64),
    Score(f32),
}

/// An exact V10 raw record. Coordinates are bin indices, oriented to the
/// user-requested axes like [`ContactRecord`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RawContactRecord {
    pub bin_x: u64,
    pub bin_y: u64,
    pub value: RawValue,
}

/// The concatenated result of a [`PreparedQuery::regions`] batch: one
/// structure-of-arrays result covering every requested region, so a caller
/// crossing an FFI or language boundary pays one call instead of one per
/// region or record.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Batch {
    /// `offsets[i]..offsets[i + 1]` is the record range for the `i`-th region.
    pub offsets: Vec<usize>,
    pub x: Vec<i64>,
    pub y: Vec<i64>,
    pub values: Vec<f32>,
}

/// A normalization vector advertised by the file footer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NormalizationEntry {
    pub normalization: Normalization,
    pub chromosome: Chromosome,
    pub unit: Unit,
    pub resolution: i32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatrixType {
    Observed,
    Oe,
    Expected,
}
impl FromStr for MatrixType {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self> {
        match s.to_ascii_lowercase().as_str() {
            "observed" => Ok(Self::Observed),
            "oe" => Ok(Self::Oe),
            "expected" => Ok(Self::Expected),
            _ => Err(Error::Argument(format!(
                "matrix type must be observed, oe, or expected (got {s})"
            ))),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Normalization(pub String);
impl Normalization {
    pub fn new(value: impl Into<String>) -> Self {
        Self(value.into().to_ascii_uppercase())
    }
    pub fn none() -> Self {
        Self("NONE".into())
    }
    pub fn is_none(&self) -> bool {
        self.0 == "NONE"
    }
    pub fn as_str(&self) -> &str {
        &self.0
    }
}
impl FromStr for Normalization {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self> {
        Ok(Self::new(s))
    }
}
impl fmt::Display for Normalization {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unit {
    BP,
    Frag,
}
impl FromStr for Unit {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self> {
        match s.to_ascii_uppercase().as_str() {
            "BP" => Ok(Self::BP),
            "FRAG" => Ok(Self::Frag),
            _ => Err(Error::Argument(format!(
                "unit must be BP or FRAG (got {s})"
            ))),
        }
    }
}
impl fmt::Display for Unit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::BP => "BP",
            Self::Frag => "FRAG",
        })
    }
}

pub struct HicFile {
    pub(crate) source: Arc<dyn RandomAccess>,
    pub path: String,
    pub version: i32,
    pub genome_id: String,
    pub attributes: BTreeMap<String, String>,
    pub chromosomes: Vec<Chromosome>,
    /// Base-pair resolutions. Kept for source compatibility; prefer
    /// [`HicFile::bp_resolutions`].
    pub resolutions: Vec<i32>,
    pub fragment_resolutions: Vec<i32>,
    master: u64,
    by_name: HashMap<String, usize>,
    normalization_entries_cache: OnceLock<Vec<NormalizationEntry>>,
    v10: Option<V10File>,
}

impl HicFile {
    pub fn open(path: impl AsRef<str>) -> Result<Self> {
        let path = path.as_ref();
        let source = open_source(path)?;
        // Read the fixed magic/version prefix in one exact range. Besides being
        // cheaper locally, this avoids a buffered reader speculatively asking
        // an HTTP server for 64 KiB once per byte while the remote length is
        // not known yet.
        let prefix = source.read_exact_at(0, 8)?;
        if &prefix[..4] != b"HIC\0" {
            return Err(Error::Invalid("magic string is not HIC".into()));
        }
        let version = i32::from_le_bytes(prefix[4..8].try_into().unwrap());
        if version < 6 {
            return Err(Error::UnsupportedVersion(version));
        }
        if version == 10 {
            let v10 = V10File::open(source.clone())?;
            let chromosomes = v10.chromosomes()?;
            let by_name = chromosomes
                .iter()
                .enumerate()
                .map(|(i, c)| (c.name.clone(), i))
                .collect();
            return Ok(Self {
                source,
                path: path.into(),
                version,
                genome_id: v10.genome(),
                attributes: v10.attributes(),
                chromosomes,
                resolutions: v10.resolutions(Unit::BP),
                fragment_resolutions: v10.resolutions(Unit::Frag),
                master: 0,
                by_name,
                normalization_entries_cache: OnceLock::new(),
                v10: Some(v10),
            });
        }
        let mut r = Reader::new(source.clone(), 8);
        let master = positive(r.i64()?, "master index")?;
        let genome_id = r.cstring()?;
        if version > 8 {
            r.i64()?;
            r.i64()?;
        }
        let attrs = nonnegative(r.i32()?, "attribute count")?;
        let mut attributes = BTreeMap::new();
        for _ in 0..attrs {
            attributes.insert(r.cstring()?, r.cstring()?);
        }
        let count = nonnegative(r.i32()?, "chromosome count")?;
        let mut chromosomes = Vec::with_capacity(count);
        let mut by_name = HashMap::with_capacity(count);
        for index in 0..count {
            let name = r.cstring()?;
            let length = if version > 8 {
                r.i64()?
            } else {
                r.i32()? as i64
            };
            by_name.insert(name.clone(), index);
            chromosomes.push(Chromosome {
                name,
                index: index as i32,
                length,
            });
        }
        let n_res = nonnegative(r.i32()?, "resolution count")?;
        let mut resolutions = Vec::with_capacity(n_res);
        for _ in 0..n_res {
            resolutions.push(r.i32()?);
        }
        let n_frag_res = nonnegative(r.i32()?, "fragment resolution count")?;
        let mut fragment_resolutions = Vec::with_capacity(n_frag_res);
        for _ in 0..n_frag_res {
            fragment_resolutions.push(r.i32()?);
        }
        Ok(Self {
            source,
            path: path.into(),
            version,
            genome_id,
            attributes,
            chromosomes,
            resolutions,
            fragment_resolutions,
            master,
            by_name,
            normalization_entries_cache: OnceLock::new(),
            v10: None,
        })
    }

    pub fn chromosome(&self, name: &str) -> Option<&Chromosome> {
        self.by_name.get(name).map(|&i| &self.chromosomes[i])
    }

    /// Genome identifier stored in the hic header, such as `hg38`.
    pub fn genome_id(&self) -> &str {
        &self.genome_id
    }

    /// On-disk hic format version.
    pub fn version(&self) -> i32 {
        self.version
    }

    /// Header attribute/value dictionary.
    pub fn attributes(&self) -> &BTreeMap<String, String> {
        &self.attributes
    }

    /// Chromosomes in their on-disk index order, including index 0 (`All`)
    /// when the file contains it.
    pub fn chromosomes(&self) -> &[Chromosome] {
        &self.chromosomes
    }

    /// Base-pair resolutions advertised in the header.
    pub fn bp_resolutions(&self) -> &[i32] {
        &self.resolutions
    }

    /// Fragment resolutions advertised in the header.
    pub fn fragment_resolutions(&self) -> &[i32] {
        &self.fragment_resolutions
    }

    /// Resolutions advertised for a particular unit.
    pub fn resolutions_for(&self, unit: Unit) -> &[i32] {
        match unit {
            Unit::BP => self.bp_resolutions(),
            Unit::Frag => self.fragment_resolutions(),
        }
    }

    /// All normalization names available anywhere in the file. `NONE` is
    /// included because raw observed data is always available without a vector.
    pub fn normalizations(&self) -> Result<Vec<Normalization>> {
        if let Some(v10) = &self.v10 {
            return Ok(v10.normalizations());
        }
        let mut names = BTreeSet::new();
        for entry in self.normalization_entries()? {
            names.insert(entry.normalization.0.clone());
        }
        Ok(normalizations_with_none(names))
    }

    /// Detailed normalization availability by chromosome, unit, and resolution.
    /// This reads only the footer index; vector payloads and contact blocks are
    /// not loaded.
    pub fn normalization_entries(&self) -> Result<&[NormalizationEntry]> {
        if self.normalization_entries_cache.get().is_none() {
            let entries = self.read_normalization_entries()?;
            let _ = self.normalization_entries_cache.set(entries);
        }
        self.normalization_entries_cache
            .get()
            .map(Vec::as_slice)
            .ok_or_else(|| Error::Invalid("normalization cache initialization failed".into()))
    }

    fn read_normalization_entries(&self) -> Result<Vec<NormalizationEntry>> {
        if let Some(v10) = &self.v10 {
            let mut out = Vec::new();
            for norm in v10.normalizations().into_iter().filter(|n| !n.is_none()) {
                for chromosome in &self.chromosomes {
                    for &unit in &[Unit::BP, Unit::Frag] {
                        for &resolution in self.resolutions_for(unit) {
                            out.push(NormalizationEntry {
                                normalization: norm.clone(),
                                chromosome: chromosome.clone(),
                                unit,
                                resolution,
                            });
                        }
                    }
                }
            }
            return Ok(out);
        }
        let mut r = Reader::new(self.source.clone(), self.master);
        if self.version > 8 {
            r.i64()?;
        } else {
            r.i32()?;
        }
        let matrices = nonnegative(r.i32()?, "matrix index count")?;
        for _ in 0..matrices {
            r.cstring()?;
            r.i64()?;
            r.i32()?;
        }
        skip_expected_maps(&mut r, self.version, false)?;
        skip_expected_maps(&mut r, self.version, true)?;
        let count = nonnegative(r.i32()?, "normalization index count")?;
        let mut entries = Vec::with_capacity(count);
        for _ in 0..count {
            let normalization = Normalization::new(r.cstring()?);
            let chromosome_index = r.i32()?;
            let unit = Unit::from_str(&r.cstring()?)?;
            let resolution = r.i32()?;
            r.i64()?;
            if self.version > 8 {
                r.i64()?;
            } else {
                r.i32()?;
            }
            let chromosome = self
                .chromosomes
                .get(usize::try_from(chromosome_index).map_err(|_| {
                    Error::Invalid("negative chromosome index in normalization entry".into())
                })?)
                .ok_or_else(|| {
                    Error::Invalid(format!(
                        "unknown chromosome index {chromosome_index} in normalization entry"
                    ))
                })?
                .clone();
            entries.push(NormalizationEntry {
                normalization,
                chromosome,
                unit,
                resolution,
            });
        }
        Ok(entries)
    }

    /// Normalizations usable for one chromosome/unit/resolution combination.
    pub fn normalizations_for(
        &self,
        chromosome: &str,
        unit: Unit,
        resolution: i32,
    ) -> Result<Vec<Normalization>> {
        let chromosome = self
            .chromosome(chromosome)
            .ok_or_else(|| Error::ChromosomeNotFound(chromosome.into()))?;
        let mut names = BTreeSet::new();
        for entry in self.normalization_entries()? {
            if entry.chromosome.index == chromosome.index
                && entry.unit == unit
                && entry.resolution == resolution
            {
                names.insert(entry.normalization.0.clone());
            }
        }
        Ok(normalizations_with_none(names))
    }

    /// The stored normalization vector for one chromosome/unit/resolution.
    /// Values are indexed by bin. `norm` must not be `NONE`.
    pub fn normalization_vector(
        &self,
        chromosome: &str,
        unit: Unit,
        resolution: i32,
        norm: &Normalization,
    ) -> Result<Vec<f64>> {
        if norm.is_none() {
            return Err(Error::Argument(
                "normalization_vector requires a normalization other than NONE".into(),
            ));
        }
        if let Some(v10) = &self.v10 {
            return v10.vector(false, chromosome, unit, resolution, norm);
        }
        let chromosome = self
            .chromosome(chromosome)
            .ok_or_else(|| Error::ChromosomeNotFound(chromosome.into()))?;
        let footer = self.read_footer(
            chromosome.index,
            chromosome.index,
            MatrixType::Observed,
            norm,
            unit,
            resolution,
        )?;
        let entry = footer
            .norm1
            .ok_or_else(|| norm_error(norm, chromosome, unit, resolution))?;
        self.read_norm(entry)
    }

    /// The stored expected-value vector for one chromosome/unit/resolution,
    /// indexed by genomic distance in bins. `norm` may be `NONE` for the raw
    /// expected vector.
    pub fn expected_vector(
        &self,
        chromosome: &str,
        unit: Unit,
        resolution: i32,
        norm: &Normalization,
    ) -> Result<Vec<f64>> {
        if let Some(v10) = &self.v10 {
            return v10.vector(true, chromosome, unit, resolution, norm);
        }
        let chromosome = self
            .chromosome(chromosome)
            .ok_or_else(|| Error::ChromosomeNotFound(chromosome.into()))?;
        let footer = self.read_footer(
            chromosome.index,
            chromosome.index,
            MatrixType::Expected,
            norm,
            unit,
            resolution,
        )?;
        if footer.expected.is_empty() {
            return Err(Error::ExpectedNotFound {
                resolution,
                unit: unit.to_string(),
            });
        }
        Ok(footer.expected)
    }

    /// Exact V10 raw records: integer counts and stored float scores are kept
    /// separate rather than routed through a lossy shared float, per
    /// [`RawValue`]. Only supported for V10 files.
    pub fn raw_records(
        &self,
        chr1: &str,
        chr2: &str,
        unit: Unit,
        resolution: i32,
    ) -> Result<Vec<RawContactRecord>> {
        let v10 = self.v10.as_ref().ok_or_else(|| {
            Error::Unsupported("exact raw records are available only for V10 files".into())
        })?;
        v10.raw_records(chr1, chr2, unit, resolution)
    }

    /// Prepare a query so repeated sub-region lookups do not reopen the
    /// footer, index, or normalization vectors for every call. See
    /// [`PreparedQuery::window`] and [`PreparedQuery::regions`].
    #[allow(clippy::too_many_arguments)] // Mirrors the established straw API.
    pub fn prepare(
        &self,
        mt: MatrixType,
        norm: Normalization,
        chr1: &str,
        chr2: &str,
        unit: Unit,
        resolution: i32,
    ) -> Result<PreparedQuery<'_>> {
        if self.v10.is_some() {
            let name = |s: &str| s.split(':').next().unwrap().to_string();
            return Ok(PreparedQuery {
                inner: PreparedInner::V10 {
                    file: self,
                    mt,
                    norm,
                    chr1: name(chr1),
                    chr2: name(chr2),
                    unit,
                    resolution,
                },
            });
        }
        let query = self.query(chr1, chr2)?;
        let zoom = self.load_zoom(query.first, query.second, mt, &norm, unit, resolution)?;
        Ok(PreparedQuery {
            inner: PreparedInner::Legacy {
                zoom,
                swapped: query.swapped,
            },
        })
    }

    pub fn records(
        &self,
        mt: MatrixType,
        norm: Normalization,
        chr1: &str,
        chr2: &str,
        unit: Unit,
        resolution: i32,
    ) -> Result<Vec<ContactRecord>> {
        if let Some(v10) = &self.v10 {
            return v10.records(mt, &norm, chr1, chr2, unit, resolution);
        }
        let query = self.query(chr1, chr2)?;
        let zoom = self.load_zoom(query.first, query.second, mt, &norm, unit, resolution)?;
        let mut records = zoom.records(query.region)?;
        if query.swapped {
            for record in &mut records {
                std::mem::swap(&mut record.bin_x, &mut record.bin_y);
            }
        }
        Ok(records)
    }

    #[allow(clippy::too_many_arguments)] // Mirrors the established straw API.
    pub fn stream_records<F>(
        &self,
        mt: MatrixType,
        norm: Normalization,
        chr1: &str,
        chr2: &str,
        unit: Unit,
        resolution: i32,
        mut callback: F,
    ) -> Result<()>
    where
        F: FnMut(ContactRecord),
    {
        if let Some(v10) = &self.v10 {
            return v10.stream_records(mt, &norm, chr1, chr2, unit, resolution, callback);
        }
        let query = self.query(chr1, chr2)?;
        let zoom = self.load_zoom(query.first, query.second, mt, &norm, unit, resolution)?;
        if query.swapped {
            let mut oriented = |mut record: ContactRecord| {
                std::mem::swap(&mut record.bin_x, &mut record.bin_y);
                callback(record);
            };
            zoom.stream(query.region, &mut oriented)
        } else {
            zoom.stream(query.region, &mut callback)
        }
    }

    pub fn matrix(
        &self,
        mt: MatrixType,
        norm: Normalization,
        chr1: &str,
        chr2: &str,
        unit: Unit,
        resolution: i32,
    ) -> Result<Vec<Vec<f32>>> {
        if self.v10.is_some() {
            let records = self.records(mt, norm, chr1, chr2, unit, resolution)?;
            let (_, ax, ay) = parse_location(chr1)?;
            let (_, bx, by) = parse_location(chr2)?;
            let rows = ((ay.unwrap_or_else(|| {
                self.chromosome(chr1.split(':').next().unwrap())
                    .unwrap()
                    .length
            }) + resolution as i64
                - 1)
                / resolution as i64
                - ax.unwrap_or(0) / resolution as i64) as usize;
            let cols = ((by.unwrap_or_else(|| {
                self.chromosome(chr2.split(':').next().unwrap())
                    .unwrap()
                    .length
            }) + resolution as i64
                - 1)
                / resolution as i64
                - bx.unwrap_or(0) / resolution as i64) as usize;
            let mut result = vec![vec![0.0; cols]; rows];
            let x0 = ax.unwrap_or(0) / resolution as i64;
            let y0 = bx.unwrap_or(0) / resolution as i64;
            for r in records {
                let x = r.bin_x as i64 / resolution as i64 - x0;
                let y = r.bin_y as i64 / resolution as i64 - y0;
                if x >= 0 && y >= 0 && (x as usize) < rows && (y as usize) < cols {
                    result[x as usize][y as usize] = r.counts;
                }
                if chr1.split(':').next() == chr2.split(':').next() {
                    let reflected_x = r.bin_y as i64 / resolution as i64 - x0;
                    let reflected_y = r.bin_x as i64 / resolution as i64 - y0;
                    if reflected_x >= 0
                        && reflected_y >= 0
                        && (reflected_x as usize) < rows
                        && (reflected_y as usize) < cols
                    {
                        result[reflected_x as usize][reflected_y as usize] = r.counts;
                    }
                }
            }
            return Ok(result);
        }
        let query = self.query(chr1, chr2)?;
        let zoom = self.load_zoom(query.first, query.second, mt, &norm, unit, resolution)?;
        let records = zoom.records(query.region)?;
        if records.is_empty() {
            return Ok(vec![vec![0.0]]);
        }
        let bins = query.region.map(|v| v / resolution as i64);
        let rows = usize::try_from(bins[1] - bins[0] + 1)
            .map_err(|_| Error::Argument("matrix row range is too large".into()))?;
        let cols = usize::try_from(bins[3] - bins[2] + 1)
            .map_err(|_| Error::Argument("matrix column range is too large".into()))?;
        let cells = rows
            .checked_mul(cols)
            .ok_or_else(|| Error::Argument("matrix dimensions overflow".into()))?;
        if cells > isize::MAX as usize / 4 {
            return Err(Error::Argument("matrix is too large to allocate".into()));
        }
        let mut matrix = vec![vec![0.0; cols]; rows];
        for rec in records {
            let row = rec.bin_x as i64 / resolution as i64 - bins[0];
            let col = rec.bin_y as i64 / resolution as i64 - bins[2];
            if row >= 0 && col >= 0 && row < rows as i64 && col < cols as i64 {
                matrix[row as usize][col as usize] = rec.counts;
            }
            if zoom.intra {
                let row = rec.bin_y as i64 / resolution as i64 - bins[0];
                let col = rec.bin_x as i64 / resolution as i64 - bins[2];
                if row >= 0 && col >= 0 && row < rows as i64 && col < cols as i64 {
                    matrix[row as usize][col as usize] = rec.counts;
                }
            }
        }
        if query.swapped {
            let mut transposed = vec![vec![0.0; rows]; cols];
            for (r, row) in matrix.into_iter().enumerate() {
                for (c, value) in row.into_iter().enumerate() {
                    transposed[c][r] = value;
                }
            }
            return Ok(transposed);
        }
        Ok(matrix)
    }

    pub fn count_records(&self, resolution: i32, inter_only: bool) -> Result<u64> {
        if let Some(v10) = &self.v10 {
            return v10.count(resolution, inter_only);
        }
        let chroms: Vec<_> = self.chromosomes.iter().filter(|c| c.index > 0).collect();
        let mut total = 0u64;
        for i in 0..chroms.len() {
            for j in (i + usize::from(inter_only))..chroms.len() {
                let z = self.load_zoom(
                    chroms[i],
                    chroms[j],
                    MatrixType::Observed,
                    &Normalization::none(),
                    Unit::BP,
                    resolution,
                )?;
                for e in z.blocks.values() {
                    total = total
                        .checked_add(record_count(&self.source, *e, self.version)?)
                        .ok_or_else(|| Error::Invalid("record count overflow".into()))?;
                }
            }
        }
        Ok(total)
    }

    pub fn chromosome_record_counts(&self, resolution: i32) -> Result<Vec<(Chromosome, u64)>> {
        if let Some(v10) = &self.v10 {
            return v10
                .chromosome_record_counts(resolution)?
                .into_iter()
                .map(|(name, count)| {
                    let chromosome = self
                        .chromosome(&name)
                        .cloned()
                        .ok_or_else(|| Error::ChromosomeNotFound(name.clone()))?;
                    Ok((chromosome, count))
                })
                .collect();
        }
        let mut result = Vec::new();
        for chromosome in self.chromosomes.iter().filter(|c| c.index > 0) {
            let zoom = self.load_zoom(
                chromosome,
                chromosome,
                MatrixType::Observed,
                &Normalization::none(),
                Unit::BP,
                resolution,
            )?;
            let mut total = 0u64;
            for entry in zoom.blocks.values() {
                total = total
                    .checked_add(record_count(&self.source, *entry, self.version)?)
                    .ok_or_else(|| Error::Invalid("record count overflow".into()))?;
            }
            result.push((chromosome.clone(), total));
        }
        Ok(result)
    }

    fn query<'a>(&'a self, a: &str, b: &str) -> Result<Query<'a>> {
        let (a_name, mut ax, mut ay) = parse_location(a)?;
        let (b_name, mut bx, mut by) = parse_location(b)?;
        let ca = self
            .chromosome(a_name)
            .ok_or_else(|| Error::ChromosomeNotFound(a_name.into()))?;
        let cb = self
            .chromosome(b_name)
            .ok_or_else(|| Error::ChromosomeNotFound(b_name.into()))?;
        if ax.is_none() {
            ax = Some(0);
            ay = Some(ca.length);
        }
        if bx.is_none() {
            bx = Some(0);
            by = Some(cb.length);
        }
        let ar = [ax.unwrap(), ay.unwrap()];
        let br = [bx.unwrap(), by.unwrap()];
        if ar[0] < 0 || br[0] < 0 || ar[1] < ar[0] || br[1] < br[0] {
            return Err(Error::Argument("invalid genomic range".into()));
        }
        if ca.index <= cb.index {
            Ok(Query {
                first: ca,
                second: cb,
                region: [ar[0], ar[1], br[0], br[1]],
                swapped: false,
            })
        } else {
            Ok(Query {
                first: cb,
                second: ca,
                region: [br[0], br[1], ar[0], ar[1]],
                swapped: true,
            })
        }
    }

    pub(crate) fn load_zoom<'a>(
        &'a self,
        c1: &Chromosome,
        c2: &Chromosome,
        mt: MatrixType,
        norm: &Normalization,
        unit: Unit,
        resolution: i32,
    ) -> Result<MatrixZoomData<'a>> {
        let footer = self.read_footer(c1.index, c2.index, mt, norm, unit, resolution)?;
        let norm1 = if norm.is_none() {
            Vec::new()
        } else {
            self.read_norm(
                footer
                    .norm1
                    .ok_or_else(|| norm_error(norm, c1, unit, resolution))?,
            )?
        };
        let norm2 = if norm.is_none() {
            Vec::new()
        } else if c1.index == c2.index {
            norm1.clone()
        } else {
            self.read_norm(
                footer
                    .norm2
                    .ok_or_else(|| norm_error(norm, c2, unit, resolution))?,
            )?
        };
        let mut r = Reader::new(self.source.clone(), footer.matrix_position);
        r.i32()?;
        r.i32()?;
        let zoom_count = nonnegative(r.i32()?, "zoom count")?;
        let mut found = None;
        for _ in 0..zoom_count {
            let zoom_unit = r.cstring()?;
            r.i32()?;
            let sum_counts = r.f32()?;
            r.f32()?;
            r.f32()?;
            r.f32()?;
            let bin_size = r.i32()?;
            let block_bin_count = r.i32()?;
            let block_column_count = r.i32()?;
            let n_blocks = nonnegative(r.i32()?, "block count")?;
            if zoom_unit == unit.to_string() && bin_size == resolution {
                let mut blocks = BTreeMap::new();
                for _ in 0..n_blocks {
                    let number = r.i32()?;
                    let position = positive(r.i64()?, "block position")?;
                    let size = positive(r.i32()? as i64, "block size")?;
                    blocks.insert(number, IndexEntry { position, size });
                }
                found = Some((sum_counts, block_bin_count, block_column_count, blocks));
                break;
            } else {
                r.skip(n_blocks as u64 * 16)?;
            }
        }
        let (sum_counts, block_bin_count, block_column_count, blocks) =
            found.ok_or_else(|| Error::ResolutionNotFound {
                resolution,
                unit: unit.to_string(),
            })?;
        let bins1 = c1.length / resolution as i64;
        let bins2 = c2.length / resolution as i64;
        let average = if c1.index == c2.index {
            0.0
        } else {
            ((sum_counts / bins1 as f32) / bins2 as f32) as f64
        };
        Ok(MatrixZoomData {
            file: self,
            intra: c1.index == c2.index,
            version: self.version,
            resolution,
            mt,
            norm: norm.clone(),
            expected: footer.expected,
            norm1,
            norm2,
            average,
            block_bin_count,
            block_column_count,
            blocks,
        })
    }

    fn read_norm(&self, entry: IndexEntry) -> Result<Vec<f64>> {
        let size = usize::try_from(entry.size)
            .map_err(|_| Error::Invalid("normalization vector is too large".into()))?;
        let data = self.source.read_exact_at(entry.position, size)?;
        let mut r = Bytes::new(&data);
        let n = if self.version > 8 {
            r.i64()?
        } else {
            r.i32()? as i64
        };
        let n = usize::try_from(n)
            .map_err(|_| Error::Invalid("invalid normalization length".into()))?;
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            out.push(if self.version > 8 {
                r.f32()? as f64
            } else {
                r.f64()?
            });
        }
        Ok(out)
    }

    fn read_footer(
        &self,
        c1: i32,
        c2: i32,
        mt: MatrixType,
        norm: &Normalization,
        unit: Unit,
        resolution: i32,
    ) -> Result<Footer> {
        let mut r = Reader::new(self.source.clone(), self.master);
        if self.version > 8 {
            r.i64()?;
        } else {
            r.i32()?;
        }
        let n = nonnegative(r.i32()?, "matrix index count")?;
        let key = format!("{c1}_{c2}");
        let mut matrix = None;
        for _ in 0..n {
            let k = r.cstring()?;
            let pos = positive(r.i64()?, "matrix position")?;
            r.i32()?;
            if k == key {
                matrix = Some(pos);
            }
        }
        let matrix_position = matrix.ok_or_else(|| Error::MatrixNotFound(key))?;
        let simple = (mt == MatrixType::Observed && norm.is_none())
            || ((mt == MatrixType::Oe || mt == MatrixType::Expected) && norm.is_none() && c1 != c2);
        if simple {
            return Ok(Footer {
                matrix_position,
                expected: Vec::new(),
                norm1: None,
                norm2: None,
            });
        }
        let mut expected = Vec::new();
        let n_exp = nonnegative(r.i32()?, "expected vector count")?;
        for _ in 0..n_exp {
            let u = r.cstring()?;
            let bs = r.i32()?;
            let nv = read_len(&mut r, self.version)?;
            let store = c1 == c2
                && mt != MatrixType::Observed
                && norm.is_none()
                && u == unit.to_string()
                && bs == resolution;
            read_expected(&mut r, self.version, nv, store, &mut expected)?;
            normalize_expected(&mut r, self.version, c1, store, &mut expected)?;
        }
        if c1 == c2 && mt != MatrixType::Observed && norm.is_none() {
            if expected.is_empty() {
                return Err(Error::ExpectedNotFound {
                    resolution,
                    unit: unit.to_string(),
                });
            }
            return Ok(Footer {
                matrix_position,
                expected,
                norm1: None,
                norm2: None,
            });
        }
        let n_exp = nonnegative(r.i32()?, "normalized expected vector count")?;
        for _ in 0..n_exp {
            let nt = r.cstring()?;
            let u = r.cstring()?;
            let bs = r.i32()?;
            let nv = read_len(&mut r, self.version)?;
            let store = c1 == c2
                && mt != MatrixType::Observed
                && nt == norm.0
                && u == unit.to_string()
                && bs == resolution;
            read_expected(&mut r, self.version, nv, store, &mut expected)?;
            normalize_expected(&mut r, self.version, c1, store, &mut expected)?;
        }
        if c1 == c2 && mt != MatrixType::Observed && !norm.is_none() && expected.is_empty() {
            return Err(Error::ExpectedNotFound {
                resolution,
                unit: unit.to_string(),
            });
        }
        let n_norm = nonnegative(r.i32()?, "normalization index count")?;
        let mut norm1 = None;
        let mut norm2 = None;
        for _ in 0..n_norm {
            let nt = r.cstring()?;
            let chr = r.i32()?;
            let u = r.cstring()?;
            let res = r.i32()?;
            let pos = positive(r.i64()?, "normalization position")?;
            let size = if self.version > 8 {
                positive(r.i64()?, "normalization size")?
            } else {
                positive(r.i32()? as i64, "normalization size")?
            };
            if nt == norm.0 && u == unit.to_string() && res == resolution {
                if chr == c1 {
                    norm1 = Some(IndexEntry {
                        position: pos,
                        size,
                    });
                }
                if chr == c2 {
                    norm2 = Some(IndexEntry {
                        position: pos,
                        size,
                    });
                }
            }
        }
        Ok(Footer {
            matrix_position,
            expected,
            norm1,
            norm2,
        })
    }
}

struct Query<'a> {
    first: &'a Chromosome,
    second: &'a Chromosome,
    region: [i64; 4],
    /// Whether `first`/`second` (and `region`) were swapped from the caller's
    /// requested order to match on-disk chromosome-index order. Callers must
    /// undo this swap on returned coordinates so results stay oriented to the
    /// axes the caller asked for, not internal storage order.
    swapped: bool,
}
struct Footer {
    matrix_position: u64,
    expected: Vec<f64>,
    norm1: Option<IndexEntry>,
    norm2: Option<IndexEntry>,
}

/// A query prepared once via [`HicFile::prepare`] and reused across repeated
/// sub-region window or batch calls, so the footer, block index, and
/// normalization vectors are not re-read for every call.
pub struct PreparedQuery<'a> {
    inner: PreparedInner<'a>,
}
enum PreparedInner<'a> {
    Legacy {
        zoom: MatrixZoomData<'a>,
        swapped: bool,
    },
    V10 {
        file: &'a HicFile,
        mt: MatrixType,
        norm: Normalization,
        chr1: String,
        chr2: String,
        unit: Unit,
        resolution: i32,
    },
}
impl PreparedQuery<'_> {
    /// A single-region window query. `region` is `[x_start, x_end, y_start,
    /// y_end]` in genomic (BP) or fragment coordinates, oriented to the axes
    /// the query was prepared with.
    pub fn window(&self, region: [i64; 4]) -> Result<Vec<ContactRecord>> {
        match &self.inner {
            PreparedInner::Legacy { zoom, swapped } => {
                let mut records = zoom.records(region)?;
                if *swapped {
                    for record in &mut records {
                        std::mem::swap(&mut record.bin_x, &mut record.bin_y);
                    }
                }
                Ok(records)
            }
            PreparedInner::V10 {
                file,
                mt,
                norm,
                chr1,
                chr2,
                unit,
                resolution,
            } => {
                let v10 = file.v10.as_ref().ok_or_else(|| {
                    Error::Invalid("prepared V10 query used against a non-V10 file".into())
                })?;
                let loc1 = format!("{chr1}:{}:{}", region[0], region[1]);
                let loc2 = format!("{chr2}:{}:{}", region[2], region[3]);
                v10.records(*mt, norm, &loc1, &loc2, *unit, *resolution)
            }
        }
    }

    /// Batched region queries as one concatenated structure-of-arrays result,
    /// so a caller pays one call per region rather than one call per record.
    pub fn regions(&self, regions: &[[i64; 4]]) -> Result<Batch> {
        let mut batch = Batch {
            offsets: Vec::with_capacity(regions.len() + 1),
            ..Batch::default()
        };
        batch.offsets.push(0);
        for region in regions {
            for record in self.window(*region)? {
                batch.x.push(record.bin_x as i64);
                batch.y.push(record.bin_y as i64);
                batch.values.push(record.counts);
            }
            batch.offsets.push(batch.x.len());
        }
        Ok(batch)
    }
}

pub(crate) struct MatrixZoomData<'a> {
    pub(crate) file: &'a HicFile,
    pub(crate) intra: bool,
    version: i32,
    resolution: i32,
    mt: MatrixType,
    norm: Normalization,
    expected: Vec<f64>,
    norm1: Vec<f64>,
    norm2: Vec<f64>,
    average: f64,
    block_bin_count: i32,
    block_column_count: i32,
    pub(crate) blocks: BTreeMap<i32, IndexEntry>,
}

impl MatrixZoomData<'_> {
    fn block_numbers(&self, region: [i64; 4]) -> Result<Vec<i32>> {
        if self.block_bin_count <= 0 || self.block_column_count <= 0 {
            return Err(Error::Invalid("invalid block geometry".into()));
        }
        let b = region.map(|v| v / self.resolution as i64);
        let mut set = BTreeSet::new();
        if self.version > 8 && self.intra {
            let low = (b[0] + b[2]) / 2 / self.block_bin_count as i64;
            let high = (b[1] + b[3]) / 2 / self.block_bin_count as i64 + 1;
            let d1 = ((1.0
                + (b[0] - b[3]).abs() as f64 / 2f64.sqrt() / self.block_bin_count as f64)
                .log2()) as i64;
            let d2 = ((1.0
                + (b[1] - b[2]).abs() as f64 / 2f64.sqrt() / self.block_bin_count as f64)
                .log2()) as i64;
            let mut near = d1.min(d2);
            if (b[0] > b[3] && b[1] < b[2]) || (b[1] > b[2] && b[0] < b[3]) {
                near = 0;
            }
            let far = d1.max(d2) + 1;
            for depth in near..=far {
                for pad in low..=high {
                    if let Ok(v) = i32::try_from(depth * self.block_column_count as i64 + pad) {
                        set.insert(v);
                    }
                }
            }
        } else {
            let c1 = b[0] / self.block_bin_count as i64;
            let c2 = (b[1] + 1) / self.block_bin_count as i64;
            let r1 = b[2] / self.block_bin_count as i64;
            let r2 = (b[3] + 1) / self.block_bin_count as i64;
            for row in r1..=r2 {
                for col in c1..=c2 {
                    if let Ok(v) = i32::try_from(row * self.block_column_count as i64 + col) {
                        set.insert(v);
                    }
                }
            }
            if self.intra {
                for row in c1..=c2 {
                    for col in r1..=r2 {
                        if let Ok(v) = i32::try_from(row * self.block_column_count as i64 + col) {
                            set.insert(v);
                        }
                    }
                }
            }
        }
        Ok(set
            .into_iter()
            .filter(|n| self.blocks.contains_key(n))
            .collect())
    }
    fn records(&self, region: [i64; 4]) -> Result<Vec<ContactRecord>> {
        let nums = self.block_numbers(region)?;
        let process = |n: &i32| -> Result<Vec<ContactRecord>> {
            let raw = read_block(&self.file.source, self.blocks[n], self.version)?;
            self.filter(raw, region)
        };
        let chunks: Result<Vec<Vec<ContactRecord>>> = nums.par_iter().map(process).collect();
        Ok(chunks?.into_iter().flatten().collect())
    }
    fn stream<F: FnMut(ContactRecord)>(&self, region: [i64; 4], callback: &mut F) -> Result<()> {
        for number in self.block_numbers(region)? {
            let raw = read_block(&self.file.source, self.blocks[&number], self.version)?;
            for record in self.filter(raw, region)? {
                callback(record);
            }
        }
        Ok(())
    }
    fn filter(&self, raw: Vec<RawRecord>, region: [i64; 4]) -> Result<Vec<ContactRecord>> {
        let mut out = Vec::new();
        for rec in raw {
            let x = rec.bin_x as i64 * self.resolution as i64;
            let y = rec.bin_y as i64 * self.resolution as i64;
            let direct = x >= region[0] && x <= region[1] && y >= region[2] && y <= region[3];
            // Cis contacts are stored above the diagonal. A below-diagonal
            // request only matches the reflected half of the stored contact,
            // so the emitted coordinates must be reflected too (below) to stay
            // oriented to the requested axes rather than storage order.
            let reflected =
                self.intra && y >= region[0] && y <= region[1] && x >= region[2] && x <= region[3];
            if !(direct || reflected) {
                continue;
            }
            let mut c = rec.counts;
            if !self.norm.is_none() {
                let a = *self
                    .norm1
                    .get(rec.bin_x as usize)
                    .ok_or_else(|| Error::Invalid("normalization vector too short".into()))?;
                let b = *self
                    .norm2
                    .get(rec.bin_y as usize)
                    .ok_or_else(|| Error::Invalid("normalization vector too short".into()))?;
                c = (c as f64 / (a * b)) as f32;
            }
            if self.mt == MatrixType::Oe {
                c = (c as f64
                    / if self.intra {
                        self.expected_value(rec.bin_x, rec.bin_y)?
                    } else {
                        self.average
                    }) as f32;
            } else if self.mt == MatrixType::Expected {
                c = (if self.intra {
                    self.expected_value(rec.bin_x, rec.bin_y)?
                } else {
                    self.average
                }) as f32;
            }
            if c.is_finite() {
                let (out_x, out_y) = if direct { (x, y) } else { (y, x) };
                out.push(ContactRecord {
                    bin_x: out_x as i32,
                    bin_y: out_y as i32,
                    counts: c,
                });
            }
        }
        Ok(out)
    }
    fn expected_value(&self, x: i32, y: i32) -> Result<f64> {
        if self.expected.is_empty() {
            return Err(Error::ExpectedNotFound {
                resolution: self.resolution,
                unit: "BP/FRAG".into(),
            });
        }
        Ok(self.expected[(x.abs_diff(y) as usize).min(self.expected.len() - 1)])
    }
}

fn parse_location(s: &str) -> Result<(&str, Option<i64>, Option<i64>)> {
    let p: Vec<_> = s.split(':').collect();
    match p.as_slice() {
        [c] => Ok((c, None, None)),
        [c, a, b] => Ok((
            c,
            Some(
                a.parse()
                    .map_err(|_| Error::Argument(format!("invalid coordinate {a}")))?,
            ),
            Some(
                b.parse()
                    .map_err(|_| Error::Argument(format!("invalid coordinate {b}")))?,
            ),
        )),
        _ => Err(Error::Argument(format!("invalid chromosome location {s}"))),
    }
}
fn positive(v: i64, name: &str) -> Result<u64> {
    u64::try_from(v).map_err(|_| Error::Invalid(format!("negative {name}")))
}
fn nonnegative(v: i32, name: &str) -> Result<usize> {
    usize::try_from(v).map_err(|_| Error::Invalid(format!("negative {name}")))
}
fn read_len(r: &mut Reader, version: i32) -> Result<u64> {
    positive(
        if version > 8 {
            r.i64()?
        } else {
            r.i32()? as i64
        },
        "vector length",
    )
}
fn skip_expected_maps(r: &mut Reader, version: i32, normalized: bool) -> Result<()> {
    let count = nonnegative(r.i32()?, "expected vector count")?;
    for _ in 0..count {
        if normalized {
            r.cstring()?;
        }
        r.cstring()?;
        r.i32()?;
        let values = read_len(r, version)?;
        let value_width = if version > 8 { 4 } else { 8 };
        r.skip(
            values
                .checked_mul(value_width)
                .ok_or_else(|| Error::Invalid("expected vector size overflow".into()))?,
        )?;
        let factors = u64::try_from(nonnegative(r.i32()?, "normalization factor count")?)
            .map_err(|_| Error::Invalid("normalization factor count overflow".into()))?;
        let factor_width = if version > 8 { 8 } else { 12 };
        r.skip(
            factors
                .checked_mul(factor_width)
                .ok_or_else(|| Error::Invalid("normalization factor size overflow".into()))?,
        )?;
    }
    Ok(())
}
fn read_expected(
    r: &mut Reader,
    version: i32,
    n: u64,
    store: bool,
    out: &mut Vec<f64>,
) -> Result<()> {
    if store {
        out.reserve(n as usize);
        for _ in 0..n {
            out.push(if version > 8 {
                r.f32()? as f64
            } else {
                r.f64()?
            });
        }
    } else {
        r.skip(n * if version > 8 { 4 } else { 8 })?;
    }
    Ok(())
}
fn normalize_expected(
    r: &mut Reader,
    version: i32,
    chr: i32,
    store: bool,
    out: &mut [f64],
) -> Result<()> {
    let n = nonnegative(r.i32()?, "normalization factor count")?;
    for _ in 0..n {
        let c = r.i32()?;
        let v = if version > 8 {
            r.f32()? as f64
        } else {
            r.f64()?
        };
        if store && c == chr {
            for x in out.iter_mut() {
                *x /= v;
            }
        }
    }
    Ok(())
}
fn norm_error(n: &Normalization, c: &Chromosome, u: Unit, r: i32) -> Error {
    Error::NormalizationNotFound {
        norm: n.0.clone(),
        chromosome: c.index,
        resolution: r,
        unit: u.to_string(),
    }
}

fn normalizations_with_none(mut names: BTreeSet<String>) -> Vec<Normalization> {
    names.remove("NONE");
    let mut result = Vec::with_capacity(names.len() + 1);
    result.push(Normalization::none());
    result.extend(names.into_iter().map(Normalization));
    result
}
