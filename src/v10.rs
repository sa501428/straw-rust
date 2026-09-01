use crate::io::RandomAccess;
use crate::{
    Chromosome, ContactRecord, Error, MatrixType, Normalization, RawContactRecord, RawValue,
    Result, Unit,
};
use std::cmp::{max, min};
use std::collections::{BTreeMap, HashSet};
use std::io::{Cursor as IoCursor, Read};
use std::sync::{Arc, Mutex};

const LIMIT: u64 = 512 * 1024 * 1024;
fn bad(s: impl Into<String>) -> Error {
    Error::Invalid(format!("V10: {}", s.into()))
}
fn req(ok: bool, s: impl Into<String>) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(bad(s))
    }
}
fn add(a: u64, b: u64) -> Result<u64> {
    a.checked_add(b).ok_or_else(|| bad("integer overflow"))
}

#[derive(Clone)]
struct Cur<'a> {
    b: &'a [u8],
    p: usize,
}
impl<'a> Cur<'a> {
    fn new(b: &'a [u8]) -> Self {
        Self { b, p: 0 }
    }
    fn left(&self) -> usize {
        self.b.len() - self.p
    }
    fn take(&mut self, n: usize) -> Result<Cur<'a>> {
        req(n <= self.left(), "truncated record")?;
        let r = Cur::new(&self.b[self.p..self.p + n]);
        self.p += n;
        Ok(r)
    }
    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?.b[0])
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.b.try_into().unwrap()))
    }
    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.take(8)?.b.try_into().unwrap()))
    }
    fn var(&mut self) -> Result<u64> {
        let mut v = 0;
        for i in 0..10 {
            let x = self.u8()?;
            req(i < 9 || x <= 1, "overflowing ULEB128")?;
            v |= u64::from(x & 127) << (i * 7);
            if x & 128 == 0 {
                req(i == 0 || x != 0, "non-canonical ULEB128")?;
                return Ok(v);
            }
        }
        Err(bad("unterminated ULEB128"))
    }
    fn zero(&mut self, n: usize) -> Result<()> {
        for _ in 0..n {
            req(self.u8()? == 0, "nonzero reserved field")?
        }
        Ok(())
    }
    fn magic(&mut self, m: &[u8; 4]) -> Result<()> {
        req(self.take(4)?.b == m, "bad record magic")
    }
    fn string(&mut self) -> Result<String> {
        let s = self.p;
        while self.p < self.b.len() && self.b[self.p] != 0 {
            self.p += 1
        }
        req(
            self.p < self.b.len() && self.p - s <= 1024 * 1024,
            "invalid string",
        )?;
        let v = std::str::from_utf8(&self.b[s..self.p])
            .map_err(|_| bad("invalid UTF-8"))?
            .into();
        self.p += 1;
        Ok(v)
    }
    fn done(&self) -> Result<()> {
        req(self.p == self.b.len(), "trailing bytes in record")
    }
}

#[derive(Clone, Copy, Default)]
struct Loc {
    pos: u64,
    len: u64,
}
fn loc(c: &mut Cur<'_>) -> Result<Loc> {
    let x = Loc {
        pos: c.u64()?,
        len: c.u64()?,
    };
    req((x.pos == 0) == (x.len == 0), "incomplete locator")?;
    Ok(x)
}
#[derive(Clone, Copy)]
struct Resolution {
    bin: u32,
    mode: u8,
    aggregation: u8,
    source: u32,
}
struct Header {
    footer: Loc,
    norm: Loc,
    expected: Loc,
    norm_expected: Loc,
    genome: String,
    attributes: Vec<(String, String)>,
    chromosomes: Vec<Chromosome>,
    resolutions: [Vec<Resolution>; 2],
    fragments: Vec<u64>,
    norms: Vec<String>,
}
impl Header {
    fn bins(&self, chr: u32, unit: u8, ri: u32) -> Result<u64> {
        let n = if unit == 1 {
            *self
                .fragments
                .get(chr as usize)
                .ok_or_else(|| bad("chromosome index"))?
        } else {
            self.chromosomes
                .get(chr as usize)
                .ok_or_else(|| bad("chromosome index"))?
                .length as u64
        };
        let b = u64::from(
            self.resolutions[unit as usize]
                .get(ri as usize)
                .ok_or_else(|| bad("resolution index"))?
                .bin,
        );
        Ok(n / b + u64::from(n % b != 0))
    }
}

fn parse_header(b: &[u8]) -> Result<Header> {
    let mut c = Cur::new(b);
    c.magic(b"HIC\0")?;
    req(c.u32()? == 10, "unsupported file version")?;
    req(
        c.u64()? == b.len() as u64 && b.len() >= 88,
        "header length mismatch",
    )?;
    let footer = loc(&mut c)?;
    req(footer.pos != 0, "missing footer")?;
    let norm = loc(&mut c)?;
    let expected = loc(&mut c)?;
    let norm_expected = loc(&mut c)?;
    c.zero(8)?;
    let genome = c.string()?;
    let n = c.u32()? as usize;
    req(n <= c.left() / 2, "attribute count")?;
    let mut attributes = Vec::with_capacity(n);
    for _ in 0..n {
        attributes.push((c.string()?, c.string()?));
    }
    let n = c.u32()? as usize;
    req(
        n > 0 && n <= c.left() / 10 && n <= i32::MAX as usize,
        "chromosome count",
    )?;
    let mut chromosomes = Vec::with_capacity(n);
    let mut names = HashSet::new();
    for i in 0..n {
        let name = c.string()?;
        let length = c.u64()?;
        req(
            !name.is_empty()
                && names.insert(name.clone())
                && length > 0
                && length <= i64::MAX as u64,
            "invalid chromosome",
        )?;
        chromosomes.push(Chromosome {
            name,
            index: i as i32,
            length: length as i64,
        });
    }
    let mut resolutions: [Vec<Resolution>; 2] = [Vec::new(), Vec::new()];
    for list in &mut resolutions {
        let n = c.u32()? as usize;
        req(n <= c.left() / 12, "resolution count")?;
        for i in 0..n {
            let r = Resolution {
                bin: c.u32()?,
                mode: c.u8()?,
                aggregation: c.u8()?,
                source: {
                    c.zero(2)?;
                    c.u32()?
                },
            };
            req(
                r.bin > 0
                    && r.bin <= i32::MAX as u32
                    && list.last().is_none_or(|x| r.bin > x.bin)
                    && r.mode <= 1
                    && r.aggregation <= 1,
                "invalid resolution",
            )?;
            if r.mode == 0 {
                req(r.source == u32::MAX, "materialized source")?
            } else {
                req(
                    r.aggregation == 1
                        && (r.source as usize) < i
                        && list[r.source as usize].mode == 0
                        && r.bin.is_multiple_of(list[r.source as usize].bin),
                    "derived source",
                )?
            }
            list.push(r)
        }
    }
    for r in &resolutions[0] {
        let sb = match r.bin {
            20 | 50 => 10,
            200 | 500 => 100,
            2000 => 1000,
            _ => 0,
        };
        if sb != 0 {
            let si = resolutions[0].iter().position(|s| s.bin == sb);
            req(
                si.is_some() && r.mode == 1 && r.source as usize == si.unwrap(),
                "mandatory BP derivation",
            )?
        }
        req(
            r.bin != 500000 || r.mode == 0,
            "500 kb must be materialized",
        )?
    }
    let mut fragments = vec![0; chromosomes.len()];
    if !resolutions[1].is_empty() {
        for (i, f) in fragments.iter_mut().enumerate() {
            let n = c.u32()? as usize;
            req(n <= c.left() / 8, "fragment count")?;
            let mut prev = 0;
            for _ in 0..n {
                let x = c.u64()?;
                req(
                    x > prev && x < chromosomes[i].length as u64,
                    "invalid fragment site",
                )?;
                prev = x
            }
            *f = n as u64 + 1
        }
    }
    let n = c.u32()? as usize;
    req(n <= c.left() / 2, "normalization count")?;
    names.clear();
    let mut norms = Vec::with_capacity(n);
    for _ in 0..n {
        let x = c.string()?;
        req(
            !x.is_empty() && x != "NONE" && names.insert(x.clone()),
            "invalid normalization",
        )?;
        norms.push(x)
    }
    c.done()?;
    let h = Header {
        footer,
        norm,
        expected,
        norm_expected,
        genome,
        attributes,
        chromosomes,
        resolutions,
        fragments,
        norms,
    };
    for u in 0..2 {
        for r in 0..h.resolutions[u].len() {
            for ch in 0..h.chromosomes.len() {
                req(
                    h.bins(ch as u32, u as u8, r as u32)? <= u32::MAX as u64,
                    "bin count exceeds uint32",
                )?
            }
        }
    }
    Ok(h)
}

#[derive(Clone)]
struct Zoom {
    unit: u8,
    mode: u8,
    value_type: u8,
    grid: u8,
    ri: u32,
    bin: u32,
    source: u32,
    b: u32,
    columns: u32,
    blocks: u32,
    sum: u64,
    occupied: u64,
    index: Loc,
}
#[derive(Clone, Copy)]
struct BlockEntry {
    number: u32,
    len: u32,
    pos: u64,
}
#[derive(Clone, Copy)]
struct Raw {
    x: u32,
    y: u32,
    value: RawValue,
}
#[derive(Default)]
struct Caches {
    zooms: BTreeMap<(u32, u32), Arc<Vec<Zoom>>>,
    indexes: BTreeMap<u64, Arc<Vec<BlockEntry>>>,
}
pub(crate) struct V10File {
    source: Arc<dyn RandomAccess>,
    header: Header,
    matrices: BTreeMap<(u32, u32), Loc>,
    caches: Mutex<Caches>,
}

fn interval(s: &Arc<dyn RandomAccess>, p: u64, n: u64) -> Result<()> {
    req(n > 0 && p.checked_add(n).is_some(), "invalid interval")?;
    if let Some(len) = s.length() {
        req(p + n <= len, "interval out of bounds")?
    }
    Ok(())
}
fn read(s: &Arc<dyn RandomAccess>, l: Loc) -> Result<Vec<u8>> {
    req(l.len <= LIMIT, "record too large")?;
    interval(s, l.pos, l.len)?;
    s.read_exact_at(
        l.pos,
        usize::try_from(l.len).map_err(|_| bad("record too large"))?,
    )
}
fn unzip(b: &[u8], n: u32) -> Result<Vec<u8>> {
    req(
        n > 0 && u64::from(n) <= LIMIT && b.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]),
        "invalid Zstandard frame",
    )?;
    let mut d =
        ruzstd::StreamingDecoder::new(IoCursor::new(b)).map_err(|e| bad(format!("zstd: {e:?}")))?;
    let mut out = Vec::with_capacity(n as usize);
    d.read_to_end(&mut out)
        .map_err(|e| bad(format!("zstd: {e}")))?;
    req(
        d.decoder.bytes_read_from_source() == b.len() as u64,
        "concatenated or trailing Zstandard data",
    )?;
    req(out.len() == n as usize, "decompressed length mismatch")?;
    Ok(out)
}

impl V10File {
    pub fn open(source: Arc<dyn RandomAccess>) -> Result<Self> {
        let p = source.read_exact_at(0, 88)?;
        let mut c = Cur::new(&p);
        c.magic(b"HIC\0")?;
        req(c.u32()? == 10, "not V10")?;
        let n = c.u64()?;
        req((88..=LIMIT).contains(&n), "header length")?;
        let header = parse_header(&source.read_exact_at(0, n as usize)?)?;
        for l in [
            header.footer,
            header.norm,
            header.expected,
            header.norm_expected,
        ] {
            if l.len > 0 {
                interval(&source, l.pos, l.len)?
            }
        }
        let b = read(&source, header.footer)?;
        let mut c = Cur::new(&b);
        c.magic(b"H10F")?;
        req(
            c.u32()? == 1 && c.u64()? == header.footer.len,
            "footer version/length",
        )?;
        let n = c.u32()? as usize;
        c.zero(4)?;
        req(n.checked_mul(24) == Some(c.left()), "footer count")?;
        let mut matrices = BTreeMap::new();
        let mut prev = None;
        for _ in 0..n {
            let k = (c.u32()?, c.u32()?);
            let l = loc(&mut c)?;
            req(
                k.0 <= k.1
                    && (k.1 as usize) < header.chromosomes.len()
                    && prev.is_none_or(|x| x < k),
                "matrix key",
            )?;
            interval(&source, l.pos, l.len)?;
            matrices.insert(k, l);
            prev = Some(k)
        }
        c.done()?;
        Ok(Self {
            source,
            header,
            matrices,
            caches: Mutex::new(Caches::default()),
        })
    }
    pub fn chromosomes(&self) -> Result<Vec<Chromosome>> {
        Ok(self.header.chromosomes.clone())
    }
    pub fn genome(&self) -> String {
        self.header.genome.clone()
    }
    pub fn resolutions(&self, u: Unit) -> Vec<i32> {
        self.header.resolutions[uid(u) as usize]
            .iter()
            .map(|r| r.bin as i32)
            .collect()
    }
    pub fn normalizations(&self) -> Vec<Normalization> {
        std::iter::once(Normalization::none())
            .chain(self.header.norms.iter().cloned().map(Normalization::new))
            .collect()
    }
    pub fn attributes(&self) -> BTreeMap<String, String> {
        self.header.attributes.iter().cloned().collect()
    }
    fn chr(&self, n: &str) -> Result<u32> {
        self.header
            .chromosomes
            .iter()
            .position(|c| c.name == n)
            .map(|x| x as u32)
            .ok_or_else(|| Error::ChromosomeNotFound(n.into()))
    }
    fn ri(&self, u: u8, b: i32) -> Result<u32> {
        self.header.resolutions[u as usize]
            .iter()
            .position(|r| r.bin == b as u32)
            .map(|x| x as u32)
            .ok_or_else(|| Error::ResolutionNotFound {
                resolution: b,
                unit: if u == 0 { "BP" } else { "FRAG" }.into(),
            })
    }
    fn ni(&self, n: &Normalization) -> Result<u32> {
        self.header
            .norms
            .iter()
            .position(|x| x == n.as_str())
            .map(|x| x as u32)
            .ok_or_else(|| bad(format!("unknown normalization {n}")))
    }
    fn zooms(&self, k: (u32, u32)) -> Result<Arc<Vec<Zoom>>> {
        if let Some(v) = self
            .caches
            .lock()
            .map_err(|_| bad("cache lock"))?
            .zooms
            .get(&k)
            .cloned()
        {
            return Ok(v);
        }
        let Some(l) = self.matrices.get(&k).copied() else {
            let v = Arc::new(vec![]);
            self.caches.lock().unwrap().zooms.insert(k, v.clone());
            return Ok(v);
        };
        let b = read(&self.source, l)?;
        let mut c = Cur::new(&b);
        c.magic(b"H10M")?;
        req(
            c.u32()? == 1 && c.u32()? == k.0 && c.u32()? == k.1,
            "matrix header",
        )?;
        let n = c.u32()? as usize;
        c.zero(4)?;
        let bp = self.header.resolutions[0].len();
        req(
            n == bp + self.header.resolutions[1].len() && n * 76 == c.left(),
            "matrix descriptors",
        )?;
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let unit = c.u8()?;
            let mode = c.u8()?;
            let aggregation = c.u8()?;
            let value_type = c.u8()?;
            let ri = c.u32()?;
            let bin = c.u32()?;
            let source = c.u32()?;
            let grid = c.u8()?;
            c.zero(3)?;
            let sum = c.u64()?;
            let occupied = c.u64()?;
            c.u32()?;
            c.u32()?;
            let block_bins = c.u32()?;
            let columns = c.u32()?;
            let index = loc(&mut c)?;
            let blocks = c.u32()?;
            c.zero(4)?;
            let eu = u8::from(i >= bp);
            let er = i - if eu == 1 { bp } else { 0 };
            let r = self.header.resolutions[eu as usize][er];
            req(
                unit == eu
                    && ri == er as u32
                    && bin == r.bin
                    && mode == r.mode
                    && aggregation == r.aggregation
                    && source == r.source,
                "descriptor mismatch",
            )?;
            req(
                value_type <= 1 && grid == u8::from(k.0 == k.1) && block_bins > 0 && columns > 0,
                "matrix geometry",
            )?;
            req(
                u64::from(columns)
                    == self
                        .header
                        .bins(k.0, unit, ri)?
                        .div_ceil(u64::from(block_bins)),
                "block columns",
            )?;
            if mode == 1 {
                req(index.pos == 0 && blocks == 0, "derived storage")?
            } else if occupied > 0 {
                req(index.pos > 0 && blocks > 0, "missing blocks")?
            } else {
                req(index.pos == 0 && blocks == 0, "empty storage")?
            }
            if index.len > 0 {
                interval(&self.source, index.pos, index.len)?
            }
            out.push(Zoom {
                unit,
                mode,
                value_type,
                grid,
                ri,
                bin,
                source,
                b: block_bins,
                columns,
                blocks,
                sum,
                occupied,
                index,
            })
        }
        c.done()?;
        for z in &out {
            if z.mode == 1 {
                let i = if z.unit == 1 { bp } else { 0 } + z.source as usize;
                req(z.value_type == out[i].value_type, "derived type")?
            }
        }
        let v = Arc::new(out);
        self.caches.lock().unwrap().zooms.insert(k, v.clone());
        Ok(v)
    }
    fn zoom(&self, k: (u32, u32), u: u8, r: u32) -> Result<Option<Zoom>> {
        let z = self.zooms(k)?;
        if z.is_empty() {
            return Ok(None);
        }
        let i = if u == 1 {
            self.header.resolutions[0].len()
        } else {
            0
        } + r as usize;
        Ok(Some(z.get(i).ok_or_else(|| bad("zoom absent"))?.clone()))
    }
    fn index(&self, z: &Zoom) -> Result<Arc<Vec<BlockEntry>>> {
        if z.index.len == 0 {
            return Ok(Arc::new(vec![]));
        }
        if let Some(v) = self
            .caches
            .lock()
            .map_err(|_| bad("cache lock"))?
            .indexes
            .get(&z.index.pos)
            .cloned()
        {
            return Ok(v);
        }
        let b = read(&self.source, z.index)?;
        let mut c = Cur::new(&b);
        c.magic(b"H10I")?;
        req(
            c.u32()? == 2 && c.u64()? == z.index.len,
            "block index header",
        )?;
        let n = c.u32()?;
        c.zero(4)?;
        req(
            n == z.blocks && z.index.len == 24 + u64::from(n) * 16,
            "block index count",
        )?;
        let mut v = Vec::with_capacity(n as usize);
        for _ in 0..n {
            let e = BlockEntry {
                number: c.u32()?,
                len: c.u32()?,
                pos: c.u64()?,
            };
            req(
                v.last().is_none_or(|p: &BlockEntry| e.number > p.number) && e.len > 16,
                "block index entry",
            )?;
            interval(&self.source, e.pos, u64::from(e.len))?;
            v.push(e)
        }
        c.done()?;
        let v = Arc::new(v);
        self.caches
            .lock()
            .unwrap()
            .indexes
            .insert(z.index.pos, v.clone());
        Ok(v)
    }

    fn block<F>(
        &self,
        b: &[u8],
        number: u32,
        z: &Zoom,
        k: (u32, u32),
        cb: &mut F,
    ) -> Result<(u64, u64)>
    where
        F: FnMut(Raw) -> Result<()>,
    {
        let mut c = Cur::new(b);
        req(c.u8()? == 1, "block version")?;
        let rep = c.u8()?;
        let mode = c.u8()?;
        let ty = c.u8()?;
        let flags = c.u8()?;
        c.zero(3)?;
        let x = c.u32()?;
        let y = c.u32()?;
        let w = c.u32()?;
        let h = c.u32()?;
        let n = c.u64()?;
        let np = c.u32()? as usize;
        let nv = c.u32()? as usize;
        req(
            rep <= 2 && mode <= 2 && ty == z.value_type && flags <= 1 && w > 0 && h > 0 && n > 0,
            "block header",
        )?;
        let cells = u64::from(w) * u64::from(h);
        let slots = if rep == 2 { cells } else { n };
        req(
            n <= cells && slots <= LIMIT / 8 && np.checked_add(nv) == Some(c.left()),
            "block streams",
        )?;
        req(
            u64::from(x) < self.header.bins(k.0, z.unit, z.ri)?
                && u64::from(y) < self.header.bins(k.1, z.unit, z.ri)?,
            "block offsets",
        )?;
        let mut ps = c.take(np)?;
        let mut vs = c.take(nv)?;
        if rep == 0 {
            req(flags == 0 && n <= np as u64, "sparse positions")?
        } else if rep == 1 || ty == 1 {
            req(flags == 1 && np as u64 == cells.div_ceil(8), "bitmap")?;
            if cells % 8 > 0 {
                req(ps.b[np - 1] >> (cells % 8) == 0, "bitmap padding")?
            }
        } else {
            req(flags == 0 && np == 0, "dense count positions")?
        }
        req(rep != 2 || mode == 2, "dense values")?;
        fn scalar(c: &mut Cur<'_>, ty: u8) -> Result<u64> {
            if ty == 1 {
                Ok(u64::from(c.u32()?))
            } else {
                c.var()
            }
        }
        let mut default = 0;
        let mut exceptions = Vec::new();
        if mode == 0 {
            default = scalar(&mut vs, ty)?
        } else if mode == 1 {
            default = scalar(&mut vs, ty)?;
            let ne = vs.var()?;
            req(
                ne > 0 && ne < slots && ne <= vs.left() as u64,
                "exception count",
            )?;
            let mut prev = 0;
            for i in 0..ne {
                let d = vs.var()?;
                req(i == 0 || d > 0, "duplicate exception")?;
                let o = if i == 0 { d } else { add(prev, d)? };
                req(o < slots, "exception range")?;
                exceptions.push(o);
                prev = o
            }
        } else {
            req(
                slots <= vs.left() as u64 / if ty == 1 { 4 } else { 1 },
                "truncated values",
            )?
        }
        let mut ei = 0;
        let mut emitted = 0;
        let mut sum = 0;
        let mut prev = 0;
        for ordinal in 0..slots {
            let value = if mode == 0 {
                default
            } else if mode == 1 && exceptions.get(ei) == Some(&ordinal) {
                ei += 1;
                let v = scalar(&mut vs, ty)?;
                req(v != default, "exception default")?;
                v
            } else if mode == 1 {
                default
            } else {
                scalar(&mut vs, ty)?
            };
            let (pos, present) = if rep == 0 {
                let d = ps.var()?;
                req(ordinal == 0 || d > 0, "duplicate cell")?;
                let p = if ordinal == 0 { d } else { add(prev, d)? };
                prev = p;
                (p, true)
            } else if rep == 1 {
                let mut p = prev;
                while p < cells && ps.b[p as usize / 8] & (1 << (p as usize % 8)) == 0 {
                    p += 1
                }
                req(p < cells, "bitmap population")?;
                prev = p + 1;
                (p, true)
            } else {
                let present = if ty == 1 {
                    ps.b[ordinal as usize / 8] & (1 << (ordinal as usize % 8)) != 0
                } else {
                    value != 0
                };
                if ty == 1 && !present {
                    req(value == 0, "absent score")?
                }
                (ordinal, present)
            };
            if !present {
                continue;
            }
            req(
                pos < cells && (rep == 2 || ty == 1 || value > 0),
                "cell/value",
            )?;
            let bx = u64::from(x) + pos % u64::from(w);
            let by = u64::from(y) + pos / u64::from(w);
            req(
                bx < self.header.bins(k.0, z.unit, z.ri)?
                    && by < self.header.bins(k.1, z.unit, z.ri)?
                    && (k.0 != k.1 || by >= bx)
                    && block_number(bx as u32, by as u32, z)? == number,
                "cell geometry",
            )?;
            let value = if ty == 1 {
                RawValue::Score(f32::from_bits(value as u32))
            } else {
                sum = add(sum, value)?;
                RawValue::Count(value)
            };
            cb(Raw {
                x: bx as u32,
                y: by as u32,
                value,
            })?;
            emitted += 1
        }
        if rep == 0 {
            ps.done()?
        } else if rep == 1 {
            req(
                ps.b.iter().map(|x| x.count_ones() as u64).sum::<u64>() == n,
                "bitmap population",
            )?
        }
        req(ei == exceptions.len(), "exception slots")?;
        vs.done()?;
        req(emitted == n, "occupied count")?;
        Ok((emitted, sum))
    }

    #[allow(clippy::too_many_arguments)]
    fn materialized<F>(
        &self,
        k: (u32, u32),
        z: &Zoom,
        x0: u64,
        x1: u64,
        y0: u64,
        y1: u64,
        mut cb: F,
    ) -> Result<()>
    where
        F: FnMut(Raw) -> Result<()>,
    {
        let index = self.index(z)?;
        let mut nb = 0;
        let mut nc = 0;
        let mut sum = 0;
        for (lo, hi) in ranges(z, x0, x1, y0, y1)? {
            let mut i = index.partition_point(|e| e.number < lo);
            while i < index.len() && index[i].number <= hi {
                let e = index[i];
                let b = self.source.read_exact_at(e.pos, e.len as usize)?;
                let mut c = Cur::new(&b);
                c.magic(b"H10B")?;
                req(c.u8()? == 1 && c.u8()? == 1, "block codec")?;
                c.zero(2)?;
                let raw = c.u32()?;
                req(c.u32()? == e.number && raw >= 40, "block record")?;
                let payload = unzip(&c.b[c.p..], raw)?;
                let (a, s) = self.block(&payload, e.number, z, k, &mut cb)?;
                nc = add(nc, a)?;
                sum = add(sum, s)?;
                nb += 1;
                i += 1
            }
        }
        if nb == index.len() as u64 {
            req(
                nb == u64::from(z.blocks)
                    && nc == z.occupied
                    && (z.value_type == 1 || sum == z.sum),
                "matrix statistics",
            )?
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn raw_canonical<F>(
        &self,
        k: (u32, u32),
        u: u8,
        ri: u32,
        x0: u64,
        x1: u64,
        y0: u64,
        y1: u64,
        mut cb: F,
    ) -> Result<()>
    where
        F: FnMut(Raw) -> Result<()>,
    {
        let Some(z) = self.zoom(k, u, ri)? else {
            return Ok(());
        };
        let inside = |x: u32, y: u32| {
            (u64::from(x) >= x0 && u64::from(x) < x1 && u64::from(y) >= y0 && u64::from(y) < y1)
                || (k.0 == k.1
                    && u64::from(y) >= x0
                    && u64::from(y) < x1
                    && u64::from(x) >= y0
                    && u64::from(x) < y1)
        };
        if z.mode == 0 {
            return self.materialized(k, &z, x0, x1, y0, y1, |r| {
                if inside(r.x, r.y) {
                    cb(r)?
                }
                Ok(())
            });
        }
        let s = self
            .zoom(k, u, z.source)?
            .ok_or_else(|| bad("derived source"))?;
        let f = u64::from(z.bin / s.bin);
        let mut counts: BTreeMap<(u32, u32), u64> = BTreeMap::new();
        let mut scores = Vec::new();
        let mut seen = HashSet::new();
        self.materialized(
            k,
            &s,
            x0 * f,
            min(x1 * f, self.header.bins(k.0, u, s.ri)?),
            y0 * f,
            min(y1 * f, self.header.bins(k.1, u, s.ri)?),
            |r| {
                let x = (u64::from(r.x) / f) as u32;
                let y = (u64::from(r.y) / f) as u32;
                if !inside(x, y) {
                    return Ok(());
                }
                req(
                    seen.insert((u64::from(r.y) << 32) | u64::from(r.x)),
                    "duplicate source cell",
                )?;
                match r.value {
                    RawValue::Count(v) => {
                        let p = counts.entry((y, x)).or_default();
                        *p = add(*p, v)?
                    }
                    RawValue::Score(v) => scores.push(((r.y, r.x), v)),
                }
                Ok(())
            },
        )?;
        if z.value_type == 1 {
            scores.sort_by_key(|x| x.0);
            let mut sums: BTreeMap<(u32, u32), f64> = BTreeMap::new();
            for ((y, x), v) in scores {
                req(v.is_finite(), "nonfinite source score")?;
                let p = sums
                    .entry(((u64::from(y) / f) as u32, (u64::from(x) / f) as u32))
                    .or_default();
                *p += f64::from(v);
                req(p.is_finite(), "score overflow")?
            }
            if x0 == 0
                && y0 == 0
                && x1 == self.header.bins(k.0, u, ri)?
                && y1 == self.header.bins(k.1, u, ri)?
            {
                req(
                    sums.len() as u64 == z.occupied,
                    "derived occupied count mismatch",
                )?;
            }
            for ((y, x), v) in sums {
                cb(Raw {
                    x,
                    y,
                    value: RawValue::Score(v as f32),
                })?
            }
        } else {
            if x0 == 0
                && y0 == 0
                && x1 == self.header.bins(k.0, u, ri)?
                && y1 == self.header.bins(k.1, u, ri)?
            {
                req(
                    counts.len() as u64 == z.occupied,
                    "derived occupied count mismatch",
                )?;
                let total = counts
                    .values()
                    .try_fold(0u64, |sum, value| add(sum, *value))?;
                req(total == z.sum, "derived sum mismatch")?;
            }
            for ((y, x), v) in counts {
                cb(Raw {
                    x,
                    y,
                    value: RawValue::Count(v),
                })?
            }
        }
        Ok(())
    }
    #[allow(clippy::too_many_arguments)]
    fn raw<F>(
        &self,
        a: u32,
        b: u32,
        u: u8,
        ri: u32,
        mut x0: u64,
        mut x1: u64,
        mut y0: u64,
        mut y1: u64,
        mut cb: F,
    ) -> Result<()>
    where
        F: FnMut(Raw) -> Result<()>,
    {
        x1 = min(x1, self.header.bins(a, u, ri)?);
        y1 = min(y1, self.header.bins(b, u, ri)?);
        if x0 >= x1 || y0 >= y1 {
            return Ok(());
        }
        let t = a > b;
        let k = if t {
            std::mem::swap(&mut x0, &mut y0);
            std::mem::swap(&mut x1, &mut y1);
            (b, a)
        } else {
            (a, b)
        };
        self.raw_canonical(k, u, ri, x0, x1, y0, y1, |mut r| {
            if t {
                std::mem::swap(&mut r.x, &mut r.y)
            }
            cb(r)
        })
    }
    fn region(&self, s: &str, u: u8, ri: u32) -> Result<Region> {
        let mut p = s.split(':');
        let chr = self.chr(p.next().unwrap_or(""))?;
        let len = if u == 1 {
            self.header.fragments[chr as usize]
        } else {
            self.header.chromosomes[chr as usize].length as u64
        };
        let (a, b) = match (p.next(), p.next(), p.next()) {
            (None, None, None) => (0, len),
            (Some(a), Some(b), None) => {
                req(
                    !a.is_empty()
                        && !b.is_empty()
                        && a.bytes().all(|x| x.is_ascii_digit())
                        && b.bytes().all(|x| x.is_ascii_digit()),
                    "invalid region",
                )?;
                let a = a.parse().map_err(|_| bad("invalid region"))?;
                let b = b.parse().map_err(|_| bad("invalid region"))?;
                req(a <= b && b <= len, "region bounds")?;
                (a, b)
            }
            _ => return Err(bad("region needs start:end")),
        };
        let bin = u64::from(self.header.resolutions[u as usize][ri as usize].bin);
        Ok(Region {
            chr,
            first: a / bin,
            last: if a == b {
                a / bin
            } else {
                b / bin + u64::from(b % bin != 0)
            },
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn vector_range(
        &self,
        kind: VKind,
        norm: u32,
        chr: u32,
        u: u8,
        ri: u32,
        begin: u64,
        wanted_end: u64,
    ) -> Result<(Vec<f64>, f64)> {
        let l = match kind {
            VKind::Norm => self.header.norm,
            VKind::Expected => self.header.expected,
            VKind::NormExpected => self.header.norm_expected,
        };
        req(l.len > 0, "vector capability absent")?;
        let b = read(&self.source, l)?;
        let mut c = Cur::new(&b);
        c.magic(match kind {
            VKind::Norm => b"NVI0",
            VKind::Expected => b"EVI0",
            VKind::NormExpected => b"NEVI",
        })?;
        req(c.u32()? == 1, "vector index version")?;
        let n = c.u32()?;
        c.zero(4)?;
        req(n as usize <= c.left() / 40, "vector entries")?;
        let mut previous = None;
        let mut result = None;
        for _ in 0..n {
            let len = c.u32()?;
            req(len >= 40, "vector entry size")?;
            let mut e = c.take(len as usize - 4)?;
            let ni = if kind == VKind::Expected { 0 } else { e.u32()? };
            let ci = if kind == VKind::Norm { e.u32()? } else { 0 };
            let eu = e.u8()?;
            e.zero(3)?;
            let eri = e.u32()?;
            let bin = e.u32()?;
            req(
                eu <= 1
                    && (eri as usize) < self.header.resolutions[eu as usize].len()
                    && bin == self.header.resolutions[eu as usize][eri as usize].bin
                    && (kind == VKind::Expected || (ni as usize) < self.header.norms.len())
                    && (kind != VKind::Norm || (ci as usize) < self.header.chromosomes.len()),
                "vector key",
            )?;
            let key = (ni, ci, eu, eri);
            req(previous.is_none_or(|x| x < key), "vector order")?;
            previous = Some(key);
            let count = e.u64()?;
            let nominal = e.u32()?;
            let chunks = e.u32()?;
            let required = if kind == VKind::Norm {
                self.header.bins(ci, eu, eri)?
            } else {
                let mut m = 0;
                for ch in 0..self.header.chromosomes.len() {
                    m = max(m, self.header.bins(ch as u32, eu, eri)?)
                }
                m
            };
            req(
                count == required
                    && (count == 0 || (nominal > 0 && chunks > 0))
                    && (count > 0 || chunks == 0),
                "vector length",
            )?;
            let matches = eu == u
                && eri == ri
                && (kind == VKind::Expected || ni == norm)
                && (kind != VKind::Norm || ci == chr);
            let mut scale = 1.0;
            if kind != VKind::Norm {
                let ns = e.u32()?;
                e.zero(4)?;
                req(ns as usize <= e.left() / 8, "scale count")?;
                let mut pc = None;
                for _ in 0..ns {
                    let ch = e.u32()?;
                    let bits = e.u32()?;
                    req(
                        (ch as usize) < self.header.chromosomes.len() && pc.is_none_or(|x| ch > x),
                        "scale key",
                    )?;
                    if ch == chr {
                        scale = f64::from(f32::from_bits(bits))
                    }
                    pc = Some(ch)
                }
            }
            req(
                u64::from(chunks) * 32 == e.left() as u64,
                "chunk descriptors",
            )?;
            let end = min(wanted_end, count);
            if matches {
                req(begin <= end && end - begin <= LIMIT / 8, "vector range")?;
                result = Some((Vec::with_capacity((end - begin) as usize), scale))
            }
            let mut next = 0;
            for _ in 0..chunks {
                let first = e.u64()?;
                let nc = e.u32()?;
                let transform = e.u8()?;
                let codec = e.u8()?;
                e.zero(2)?;
                let pos = e.u64()?;
                let stored = e.u32()?;
                let raw = e.u32()?;
                req(
                    first == next
                        && nc > 0
                        && u64::from(nc) * 4 == u64::from(raw)
                        && transform <= 2
                        && codec == 1
                        && stored > 16,
                    "chunk descriptor",
                )?;
                next = add(next, u64::from(nc))?;
                req(next <= count, "chunk range")?;
                interval(&self.source, pos, u64::from(stored))?;
                if !matches || next <= begin || first >= end {
                    continue;
                }
                let cb = self.source.read_exact_at(pos, stored as usize)?;
                let mut ch = Cur::new(&cb);
                ch.magic(b"H10V")?;
                req(ch.u8()? == codec && ch.u8()? == transform, "chunk codec")?;
                ch.zero(2)?;
                req(ch.u32()? == raw && ch.u32()? == nc, "chunk size")?;
                let data = unzip(&ch.b[ch.p..], raw)?;
                let mut words = Cur::new(&data);
                let mut prev = 0;
                for i in 0..nc {
                    let mut bits = if transform == 1 {
                        let mut v = 0;
                        for lane in 0..4 {
                            v |= u32::from(data[lane * nc as usize + i as usize]) << (lane * 8)
                        }
                        v
                    } else {
                        words.u32()?
                    };
                    if transform == 2 && i > 0 {
                        bits ^= prev
                    }
                    prev = bits;
                    let at = first + u64::from(i);
                    if at >= begin && at < end {
                        result
                            .as_mut()
                            .unwrap()
                            .0
                            .push(f64::from(f32::from_bits(bits)))
                    }
                }
            }
            req(next == count, "vector coverage")?;
            e.done()?
        }
        c.done()?;
        result.ok_or_else(|| bad("vector capability absent"))
    }

    pub fn vector(
        &self,
        expected: bool,
        chromosome: &str,
        unit: Unit,
        resolution: i32,
        norm: &Normalization,
    ) -> Result<Vec<f64>> {
        let ch = self.chr(chromosome)?;
        let u = uid(unit);
        let ri = self.ri(u, resolution)?;
        if !expected && norm.is_none() {
            return Ok(vec![1.0; self.header.bins(ch, u, ri)? as usize]);
        }
        let kind = if expected {
            if norm.is_none() {
                VKind::Expected
            } else {
                VKind::NormExpected
            }
        } else {
            VKind::Norm
        };
        let ni = if norm.is_none() { 0 } else { self.ni(norm)? };
        let (mut v, scale) = self.vector_range(kind, ni, ch, u, ri, 0, u64::MAX)?;
        if expected {
            req(scale != 0.0 && scale.is_finite(), "expected scale")?;
            for x in &mut v {
                *x /= scale
            }
        }
        Ok(v)
    }
    pub fn raw_records(
        &self,
        a: &str,
        b: &str,
        unit: Unit,
        resolution: i32,
    ) -> Result<Vec<RawContactRecord>> {
        let u = uid(unit);
        let ri = self.ri(u, resolution)?;
        let x = self.region(a, u, ri)?;
        let y = self.region(b, u, ri)?;
        let mut out = Vec::new();
        self.raw(
            x.chr,
            y.chr,
            u,
            ri,
            x.first,
            x.last,
            y.first,
            y.last,
            |mut r| {
                if x.chr == y.chr
                    && !(u64::from(r.x) >= x.first
                        && u64::from(r.x) < x.last
                        && u64::from(r.y) >= y.first
                        && u64::from(r.y) < y.last)
                {
                    std::mem::swap(&mut r.x, &mut r.y)
                }
                out.push(RawContactRecord {
                    bin_x: u64::from(r.x),
                    bin_y: u64::from(r.y),
                    value: r.value,
                });
                Ok(())
            },
        )?;
        Ok(out)
    }
    #[allow(clippy::too_many_arguments)]
    fn stream<F>(
        &self,
        mt: MatrixType,
        norm: &Normalization,
        a: &str,
        b: &str,
        unit: Unit,
        resolution: i32,
        mut emit: F,
    ) -> Result<()>
    where
        F: FnMut(ContactRecord),
    {
        let u = uid(unit);
        let ri = self.ri(u, resolution)?;
        let x = self.region(a, u, ri)?;
        let y = self.region(b, u, ri)?;
        if x.first == x.last || y.first == y.last {
            return Ok(());
        }
        let ni = if norm.is_none() { 0 } else { self.ni(norm)? };
        let mut nx = Vec::new();
        let mut ny = Vec::new();
        if !norm.is_none() && mt != MatrixType::Expected {
            nx = self
                .vector_range(VKind::Norm, ni, x.chr, u, ri, x.first, x.last)?
                .0;
            ny = self
                .vector_range(VKind::Norm, ni, y.chr, u, ri, y.first, y.last)?
                .0
        }
        let mut ev = Vec::new();
        let mut scale = 1.0;
        let mut eb = 0;
        if mt != MatrixType::Observed {
            req(x.chr == y.chr, "expected only for cis")?;
            eb = if x.last <= y.first {
                y.first - x.last + 1
            } else if y.last <= x.first {
                x.first - y.last + 1
            } else {
                0
            };
            let ee = max(
                if x.last > y.first {
                    x.last - y.first
                } else {
                    y.first - x.last + 2
                },
                if y.last > x.first {
                    y.last - x.first
                } else {
                    x.first - y.last + 2
                },
            );
            let pair = self.vector_range(
                if norm.is_none() {
                    VKind::Expected
                } else {
                    VKind::NormExpected
                },
                ni,
                x.chr,
                u,
                ri,
                eb,
                ee,
            )?;
            ev = pair.0;
            scale = pair.1
        }
        self.raw(
            x.chr,
            y.chr,
            u,
            ri,
            x.first,
            x.last,
            y.first,
            y.last,
            |mut r| {
                if x.chr == y.chr
                    && !(u64::from(r.x) >= x.first
                        && u64::from(r.x) < x.last
                        && u64::from(r.y) >= y.first
                        && u64::from(r.y) < y.last)
                {
                    std::mem::swap(&mut r.x, &mut r.y)
                }
                let mut value = match r.value {
                    RawValue::Count(v) => v as f64,
                    RawValue::Score(v) => f64::from(v),
                };
                if !norm.is_none() && mt != MatrixType::Expected {
                    let a = nx[(u64::from(r.x) - x.first) as usize];
                    let b = ny[(u64::from(r.y) - y.first) as usize];
                    if a == 0.0 || b == 0.0 || !a.is_finite() || !b.is_finite() {
                        return Ok(());
                    }
                    value /= a * b
                }
                if mt != MatrixType::Observed {
                    let d = u64::from(r.x.abs_diff(r.y));
                    if d < eb || d - eb >= ev.len() as u64 || scale == 0.0 || !scale.is_finite() {
                        return Ok(());
                    }
                    let e = ev[(d - eb) as usize] / scale;
                    if e == 0.0 || !e.is_finite() {
                        return Ok(());
                    }
                    value = if mt == MatrixType::Oe { value / e } else { e }
                }
                let bx = u64::from(r.x) * resolution as u64;
                let by = u64::from(r.y) * resolution as u64;
                req(
                    bx <= i32::MAX as u64 && by <= i32::MAX as u64,
                    "contact coordinate overflow",
                )?;
                emit(ContactRecord {
                    bin_x: bx as i32,
                    bin_y: by as i32,
                    counts: value as f32,
                });
                Ok(())
            },
        )?;
        Ok(())
    }
    pub fn records(
        &self,
        mt: MatrixType,
        norm: &Normalization,
        a: &str,
        b: &str,
        unit: Unit,
        resolution: i32,
    ) -> Result<Vec<ContactRecord>> {
        let mut records = Vec::new();
        self.stream(mt, norm, a, b, unit, resolution, |record| {
            records.push(record)
        })?;
        Ok(records)
    }
    #[allow(clippy::too_many_arguments)]
    pub fn stream_records<F>(
        &self,
        mt: MatrixType,
        norm: &Normalization,
        a: &str,
        b: &str,
        unit: Unit,
        resolution: i32,
        emit: F,
    ) -> Result<()>
    where
        F: FnMut(ContactRecord),
    {
        self.stream(mt, norm, a, b, unit, resolution, emit)
    }
    pub fn count(&self, resolution: i32, inter: bool) -> Result<u64> {
        let ri = self.ri(0, resolution)?;
        let mut total = 0;
        for a in 0..self.header.chromosomes.len() as u32 {
            if self.header.chromosomes[a as usize]
                .name
                .eq_ignore_ascii_case("all")
            {
                continue;
            }
            for b in a..self.header.chromosomes.len() as u32 {
                if inter && a == b {
                    continue;
                }
                if let Some(z) = self.zoom((a, b), 0, ri)? {
                    total = add(total, z.occupied)?
                }
            }
        }
        Ok(total)
    }
    pub fn chromosome_record_counts(&self, resolution: i32) -> Result<Vec<(String, u64)>> {
        let ri = self.ri(0, resolution)?;
        let mut out = Vec::new();
        for ch in 0..self.header.chromosomes.len() as u32 {
            let name = &self.header.chromosomes[ch as usize].name;
            if name.eq_ignore_ascii_case("all") {
                continue;
            }
            out.push((
                name.clone(),
                self.zoom((ch, ch), 0, ri)?.map_or(0, |z| z.occupied),
            ))
        }
        Ok(out)
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum VKind {
    Norm,
    Expected,
    NormExpected,
}
struct Region {
    chr: u32,
    first: u64,
    last: u64,
}
fn uid(u: Unit) -> u8 {
    u8::from(u == Unit::Frag)
}
fn depth(d: u64, b: u32) -> u32 {
    let lhs = u128::from(d) * u128::from(d);
    let scale = 2 * u128::from(b) * u128::from(b);
    let mut x = 0;
    while x < 32 {
        let t = (1u128 << (x + 1)) - 1;
        if t * t > lhs / scale {
            break;
        }
        x += 1
    }
    x
}
fn block_number(x: u32, y: u32, z: &Zoom) -> Result<u32> {
    let n = if z.grid == 0 {
        u64::from(y / z.b) * u64::from(z.columns) + u64::from(x / z.b)
    } else {
        u64::from(depth(u64::from(y) - u64::from(x), z.b)) * u64::from(z.columns)
            + (u64::from(x) + u64::from(y)) / (2 * u64::from(z.b))
    };
    u32::try_from(n).map_err(|_| bad("block number"))
}
fn ranges(z: &Zoom, x0: u64, x1: u64, y0: u64, y1: u64) -> Result<Vec<(u32, u32)>> {
    let mut out = Vec::new();
    if x0 >= x1 || y0 >= y1 {
        return Ok(out);
    }
    let mut push = |a: u64, b: u64| {
        if a <= u32::MAX as u64 {
            out.push((a as u32, min(b, u32::MAX as u64) as u32))
        }
    };
    if z.grid == 0 {
        let fc = x0 / u64::from(z.b);
        let lc = (x1 - 1) / u64::from(z.b);
        for row in y0 / u64::from(z.b)..=(y1 - 1) / u64::from(z.b) {
            push(
                row * u64::from(z.columns) + fc,
                row * u64::from(z.columns) + lc,
            )
        }
    } else {
        let lx = x1 - 1;
        let ly = y1 - 1;
        let nearest = if x1 <= y0 {
            y0 - lx
        } else if y1 <= x0 {
            x0 - ly
        } else {
            0
        };
        let far = max(ly.abs_diff(x0), lx.abs_diff(y0));
        let d = 2 * u64::from(z.b);
        for dep in depth(nearest, z.b)..=depth(far, z.b) {
            push(
                u64::from(dep) * u64::from(z.columns) + (x0 + y0) / d,
                u64::from(dep) * u64::from(z.columns) + (lx + ly) / d,
            )
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn leb() {
        let mut c = Cur::new(&[0xac, 2]);
        assert_eq!(c.var().unwrap(), 300);
        assert!(Cur::new(&[0x80, 0]).var().is_err())
    }
}
