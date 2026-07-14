use crate::{Error, Result};
use reqwest::blocking::Client;
use reqwest::header::RANGE;
use std::fs::File;
use std::sync::Arc;

pub(crate) trait RandomAccess: Send + Sync {
    fn read_exact_at(&self, offset: u64, len: usize) -> Result<Vec<u8>>;
    fn length(&self) -> Option<u64> {
        None
    }
}

pub(crate) fn open_source(path: &str) -> Result<Arc<dyn RandomAccess>> {
    if path.starts_with("http://") || path.starts_with("https://") {
        Ok(Arc::new(HttpSource {
            client: Client::builder().user_agent("straw-rust").build()?,
            url: path.to_owned(),
        }))
    } else {
        let file = File::open(path)?;
        let length = file.metadata()?.len();
        Ok(Arc::new(LocalSource { file, length }))
    }
}

struct LocalSource {
    file: File,
    length: u64,
}

impl RandomAccess for LocalSource {
    fn read_exact_at(&self, offset: u64, len: usize) -> Result<Vec<u8>> {
        let mut buf = vec![0; len];
        #[cfg(unix)]
        {
            use std::os::unix::fs::FileExt;
            let mut done = 0;
            while done < len {
                let n = self.file.read_at(&mut buf[done..], offset + done as u64)?;
                if n == 0 {
                    return Err(std::io::Error::from(std::io::ErrorKind::UnexpectedEof).into());
                }
                done += n;
            }
        }
        #[cfg(not(unix))]
        {
            use std::io::{Read, Seek, SeekFrom};
            let mut f = self.file.try_clone()?;
            f.seek(SeekFrom::Start(offset))?;
            f.read_exact(&mut buf)?;
        }
        Ok(buf)
    }
    fn length(&self) -> Option<u64> {
        Some(self.length)
    }
}

struct HttpSource {
    client: Client,
    url: String,
}

impl RandomAccess for HttpSource {
    fn read_exact_at(&self, offset: u64, len: usize) -> Result<Vec<u8>> {
        if len == 0 {
            return Ok(Vec::new());
        }
        let end = offset
            .checked_add(len as u64 - 1)
            .ok_or_else(|| Error::Invalid("range overflow".into()))?;
        let response = self
            .client
            .get(&self.url)
            .header(RANGE, format!("bytes={offset}-{end}"))
            .send()?
            .error_for_status()?;
        if response.status() != reqwest::StatusCode::PARTIAL_CONTENT {
            return Err(Error::Invalid(
                "HTTP server does not support byte ranges".into(),
            ));
        }
        let bytes = response.bytes()?;
        if bytes.len() < len {
            return Err(Error::Invalid(format!(
                "short HTTP range: wanted {len}, got {}",
                bytes.len()
            )));
        }
        Ok(bytes[..len].to_vec())
    }
}

/// Buffered cursor over a random-access source. Sequential metadata parsing
/// normally costs one range read per 64 KiB rather than one per scalar.
pub(crate) struct Reader {
    source: Arc<dyn RandomAccess>,
    pos: u64,
    start: u64,
    cache: Vec<u8>,
}

impl Reader {
    pub fn new(source: Arc<dyn RandomAccess>, pos: u64) -> Self {
        Self {
            source,
            pos,
            start: 0,
            cache: Vec::new(),
        }
    }
    pub fn skip(&mut self, n: u64) -> Result<()> {
        self.pos = self
            .pos
            .checked_add(n)
            .ok_or_else(|| Error::Invalid("offset overflow".into()))?;
        Ok(())
    }
    fn take(&mut self, n: usize) -> Result<&[u8]> {
        let cached_end = self.start + self.cache.len() as u64;
        if self.pos < self.start || self.pos + n as u64 > cached_end {
            self.start = self.pos;
            let desired = n.max(64 * 1024);
            let desired = self.source.length().map_or(desired, |len| {
                desired.min(len.saturating_sub(self.pos) as usize)
            });
            if desired < n {
                return Err(std::io::Error::from(std::io::ErrorKind::UnexpectedEof).into());
            }
            self.cache = self
                .source
                .read_exact_at(self.pos, desired)
                .or_else(|_| self.source.read_exact_at(self.pos, n))?;
        }
        let i = (self.pos - self.start) as usize;
        self.pos += n as u64;
        Ok(&self.cache[i..i + n])
    }
    pub fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }
    pub fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    pub fn i64(&mut self) -> Result<i64> {
        Ok(i64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    pub fn f32(&mut self) -> Result<f32> {
        Ok(f32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    pub fn f64(&mut self) -> Result<f64> {
        Ok(f64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    pub fn cstring(&mut self) -> Result<String> {
        let mut bytes = Vec::new();
        loop {
            let b = self.u8()?;
            if b == 0 {
                break;
            }
            bytes.push(b);
        }
        String::from_utf8(bytes).map_err(|e| Error::Invalid(format!("non-UTF8 string: {e}")))
    }
}

pub(crate) struct Bytes<'a> {
    data: &'a [u8],
    pos: usize,
}
impl<'a> Bytes<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self
            .pos
            .checked_add(n)
            .ok_or_else(|| Error::CorruptBlock("offset overflow".into()))?;
        let v = self
            .data
            .get(self.pos..end)
            .ok_or_else(|| Error::CorruptBlock("unexpected end".into()))?;
        self.pos = end;
        Ok(v)
    }
    pub fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }
    pub fn i16(&mut self) -> Result<i16> {
        Ok(i16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }
    pub fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    pub fn i64(&mut self) -> Result<i64> {
        Ok(i64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    pub fn f32(&mut self) -> Result<f32> {
        Ok(f32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    pub fn f64(&mut self) -> Result<f64> {
        Ok(f64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
}
