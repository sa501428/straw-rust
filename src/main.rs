use std::env;
use std::io::{self, BufWriter, Write};
use std::process::ExitCode;
use std::str::FromStr;
use straw::{ContactFilter, DumpOptions, MatrixType, Normalization, Unit};

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("Error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> straw::Result<()> {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.first().map(String::as_str) == Some("dump") {
        return run_dump(&args[1..]);
    }
    if args.len() != 6 && args.len() != 7 {
        return Err(straw::Error::Argument("Usage: straw [observed/oe/expected] <NONE/VC/VC_SQRT/KR> <hicFile> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG/MATRIX> <binsize>".into()));
    }
    let (mt, o) = if args.len() == 7 {
        (MatrixType::from_str(&args[0])?, 1)
    } else {
        (MatrixType::Observed, 0)
    };
    let norm = Normalization::new(&args[o]);
    let file = &args[o + 1];
    let c1 = &args[o + 2];
    let c2 = &args[o + 3];
    let unit_arg = args[o + 4].to_ascii_uppercase();
    let resolution: i32 = args[o + 5]
        .parse()
        .map_err(|_| straw::Error::Argument("invalid binsize".into()))?;
    if resolution <= 0 {
        return Err(straw::Error::Argument("binsize must be positive".into()));
    }
    let mut out = BufWriter::with_capacity(1024 * 1024, io::stdout().lock());
    if unit_arg == "MATRIX" {
        for row in straw::straw_as_matrix(mt, norm, file, c1, c2, Unit::BP, resolution)? {
            for value in row {
                write!(out, "{}\t", format_g(value))?;
            }
            writeln!(out)?;
        }
    } else {
        let unit = Unit::from_str(&unit_arg)?;
        straw::straw_stream(mt, norm, file, c1, c2, unit, resolution, |r| {
            let _ = writeln!(out, "{}\t{}\t{}", r.bin_x, r.bin_y, format_g(r.counts));
        })?;
    }
    out.flush()?;
    Ok(())
}

fn run_dump(a: &[String]) -> straw::Result<()> {
    if a.len() != 7 && a.len() != 8 {
        return Err(straw::Error::Argument("Usage: straw dump <observed/oe/expected> <NONE/VC/VC_SQRT/KR> <hicFile> <BP/FRAG> <binsize> <outputFile> <compressed> [-intra-short|-intra-long|-inter|-intra]".into()));
    }
    let compressed = matches!(
        a[6].to_ascii_lowercase().as_str(),
        "1" | "true" | "compressed"
    );
    let filter = if let Some(v) = a.get(7) {
        let v = v.to_ascii_lowercase();
        if v.contains("inter") {
            ContactFilter::Inter
        } else if v.contains("intra") && v.contains("short") {
            ContactFilter::IntraShort
        } else if v.contains("intra") && v.contains("long") {
            ContactFilter::IntraLong
        } else if v.contains("intra") {
            ContactFilter::Intra
        } else {
            ContactFilter::All
        }
    } else {
        ContactFilter::All
    };
    let opts = DumpOptions {
        matrix_type: MatrixType::from_str(&a[0])?,
        normalization: Normalization::new(&a[1]),
        unit: Unit::from_str(&a[3])?,
        resolution: a[4]
            .parse()
            .map_err(|_| straw::Error::Argument("invalid binsize".into()))?,
        compressed,
        filter,
    };
    straw::dump(&a[2], &a[5], &opts)
}

fn format_g(value: f32) -> String {
    let v = value as f64;
    if v == 0.0 {
        return "0".into();
    }
    let av = v.abs();
    let exp = av.log10().floor() as i32;
    if !(-4..14).contains(&exp) {
        let s = format!("{:.13e}", v);
        let (m, e) = s.split_once('e').unwrap();
        let m = m.trim_end_matches('0').trim_end_matches('.');
        let ei: i32 = e.parse().unwrap();
        format!("{m}e{ei:+03}")
    } else {
        let decimals = (13 - exp).max(0) as usize;
        format!("{v:.decimals$}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string()
    }
}
