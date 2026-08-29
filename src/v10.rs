use crate::{Chromosome, ContactRecord, Error, MatrixType, Normalization, Result, Unit};
use std::ffi::{c_char, c_int, c_void, CStr, CString};
use std::sync::Mutex;

#[repr(C)]
struct CRecord {
    x: i32,
    y: i32,
    value: f32,
}
type Callback = unsafe extern "C" fn(*mut c_void, CRecord);
extern "C" {
    fn straw_v10_error() -> *const c_char;
    fn straw_v10_open(path: *const c_char) -> *mut c_void;
    fn straw_v10_close(file: *mut c_void);
    fn straw_v10_chromosome_count(file: *mut c_void) -> usize;
    fn straw_v10_chromosome(
        file: *mut c_void,
        i: usize,
        name: *mut *const c_char,
        index: *mut i32,
        length: *mut i64,
    ) -> c_int;
    fn straw_v10_string_free(value: *const c_char);
    fn straw_v10_genome(file: *mut c_void) -> *const c_char;
    fn straw_v10_resolution_count(file: *mut c_void, frag: c_int) -> usize;
    fn straw_v10_resolution(file: *mut c_void, frag: c_int, i: usize) -> i32;
    fn straw_v10_norm_count(file: *mut c_void) -> usize;
    fn straw_v10_norm(file: *mut c_void, i: usize) -> *const c_char;
    fn straw_v10_attribute_count(file: *mut c_void) -> usize;
    fn straw_v10_attribute(file: *mut c_void, i: usize, value: c_int) -> *const c_char;
    fn straw_v10_stream(
        file: *mut c_void,
        mt: *const c_char,
        norm: *const c_char,
        a: *const c_char,
        b: *const c_char,
        unit: *const c_char,
        resolution: i32,
        context: *mut c_void,
        callback: Callback,
    ) -> c_int;
    fn straw_v10_count(file: *mut c_void, resolution: i32, inter: c_int) -> i64;
}

pub(crate) struct V10File {
    raw: usize,
    lock: Mutex<()>,
}
impl Drop for V10File {
    fn drop(&mut self) {
        unsafe { straw_v10_close(self.raw as *mut c_void) }
    }
}
impl V10File {
    pub fn open(path: &str) -> Result<Self> {
        let path = CString::new(path).map_err(|_| Error::Argument("path contains NUL".into()))?;
        let raw = unsafe { straw_v10_open(path.as_ptr()) };
        if raw.is_null() {
            Err(last_error())
        } else {
            Ok(Self {
                raw: raw as usize,
                lock: Mutex::new(()),
            })
        }
    }
    pub fn chromosomes(&self) -> Result<Vec<Chromosome>> {
        let _guard = self
            .lock
            .lock()
            .map_err(|_| Error::Invalid("V10 reader lock poisoned".into()))?;
        let n = unsafe { straw_v10_chromosome_count(self.raw as *mut c_void) };
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let mut name = std::ptr::null();
            let mut index = 0;
            let mut length = 0;
            if unsafe {
                straw_v10_chromosome(
                    self.raw as *mut c_void,
                    i,
                    &mut name,
                    &mut index,
                    &mut length,
                )
            } == 0
            {
                return Err(last_error());
            }
            let value = unsafe { CStr::from_ptr(name) }
                .to_string_lossy()
                .into_owned();
            unsafe { straw_v10_string_free(name) };
            out.push(Chromosome {
                name: value,
                index,
                length,
            });
        }
        Ok(out)
    }
    pub fn genome(&self) -> String {
        let _guard = self.lock.lock().unwrap();
        unsafe { CStr::from_ptr(straw_v10_genome(self.raw as *mut c_void)) }
            .to_string_lossy()
            .into_owned()
    }
    pub fn resolutions(&self, unit: Unit) -> Vec<i32> {
        let _guard = self.lock.lock().unwrap();
        let frag = i32::from(unit == Unit::Frag);
        let n = unsafe { straw_v10_resolution_count(self.raw as *mut c_void, frag) };
        (0..n)
            .map(|i| unsafe { straw_v10_resolution(self.raw as *mut c_void, frag, i) })
            .collect()
    }
    pub fn normalizations(&self) -> Vec<Normalization> {
        let _guard = self.lock.lock().unwrap();
        let n = unsafe { straw_v10_norm_count(self.raw as *mut c_void) };
        let mut out = vec![Normalization::none()];
        for i in 0..n {
            out.push(Normalization::new(
                unsafe { CStr::from_ptr(straw_v10_norm(self.raw as *mut c_void, i)) }
                    .to_string_lossy(),
            ));
        }
        out
    }
    pub fn attributes(&self) -> std::collections::BTreeMap<String, String> {
        let _guard = self.lock.lock().unwrap();
        let n = unsafe { straw_v10_attribute_count(self.raw as *mut c_void) };
        (0..n)
            .map(|i| {
                let key =
                    unsafe { CStr::from_ptr(straw_v10_attribute(self.raw as *mut c_void, i, 0)) }
                        .to_string_lossy()
                        .into_owned();
                let value =
                    unsafe { CStr::from_ptr(straw_v10_attribute(self.raw as *mut c_void, i, 1)) }
                        .to_string_lossy()
                        .into_owned();
                (key, value)
            })
            .collect()
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
        unsafe extern "C" fn collect(context: *mut c_void, r: CRecord) {
            (&mut *(context as *mut Vec<ContactRecord>)).push(ContactRecord {
                bin_x: r.x,
                bin_y: r.y,
                counts: r.value,
            });
        }
        let strings = [matrix_name(mt), norm.as_str(), a, b, unit_name(unit)]
            .map(|s| CString::new(s).map_err(|_| Error::Argument("argument contains NUL".into())))
            .into_iter()
            .collect::<Result<Vec<_>>>()?;
        let mut out = Vec::new();
        let _guard = self
            .lock
            .lock()
            .map_err(|_| Error::Invalid("V10 reader lock poisoned".into()))?;
        let ok = unsafe {
            straw_v10_stream(
                self.raw as *mut c_void,
                strings[0].as_ptr(),
                strings[1].as_ptr(),
                strings[2].as_ptr(),
                strings[3].as_ptr(),
                strings[4].as_ptr(),
                resolution,
                &mut out as *mut _ as *mut c_void,
                collect,
            )
        };
        if ok == 0 {
            Err(last_error())
        } else {
            Ok(out)
        }
    }
    pub fn count(&self, resolution: i32, inter: bool) -> Result<u64> {
        let _guard = self
            .lock
            .lock()
            .map_err(|_| Error::Invalid("V10 reader lock poisoned".into()))?;
        let n = unsafe { straw_v10_count(self.raw as *mut c_void, resolution, i32::from(inter)) };
        if n < 0 {
            Err(last_error())
        } else {
            Ok(n as u64)
        }
    }
}
fn matrix_name(mt: MatrixType) -> &'static str {
    match mt {
        MatrixType::Observed => "observed",
        MatrixType::Oe => "oe",
        MatrixType::Expected => "expected",
    }
}
fn unit_name(unit: Unit) -> &'static str {
    match unit {
        Unit::BP => "BP",
        Unit::Frag => "FRAG",
    }
}
fn last_error() -> Error {
    Error::Invalid(
        unsafe { CStr::from_ptr(straw_v10_error()) }
            .to_string_lossy()
            .into_owned(),
    )
}
