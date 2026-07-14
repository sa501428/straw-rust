# straw-rust

A native Rust port of the C++ straw reader for indexed access to large
Juicebox `.hic` files. It provides both a Rust library and a command-line tool.

## Features

- `.hic` versions 6 and newer, including v9 intra-chromosomal block geometry
  and v10 delta/all-one block encodings
- zlib and zstd block compression
- local files and HTTP(S) byte-range sources
- `observed`, `oe`, and `expected` matrices
- arbitrary normalization names (`NONE`, `VC`, `VC_SQRT`, `KR`, `SCALE`, etc.)
- BP and fragment units, sparse streaming, sparse vectors, and dense matrices
- C++-compatible compressed or uncompressed `HICSLICE` dumps and filters
- genome-wide and per-chromosome record counting

Local data is read with positional I/O: the complete `.hic` file is never
loaded into memory, concurrent block reads do not share a seek cursor, and
independent blocks are decompressed with Rayon. The streaming API retains only
one decompressed block at a time.

## Build

```bash
cargo build --release
```

The executable is `target/release/straw`.

## CLI

The CLI accepts the same forms as the bundled C++ implementation:

```text
straw [observed/oe/expected] <NONE/VC/VC_SQRT/KR> <hicFile> \
      <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG/MATRIX> <binsize>

straw dump <observed/oe/expected> <normalization> <hicFile> <BP/FRAG> \
      <binsize> <outputFile> <compressed> \
      [-intra-short|-intra-long|-inter|-intra]
```

The matrix type is optional and defaults to `observed`, as in C++.

```bash
target/release/straw observed NONE sample.hic \
  chr1:0:1000000 chr1:0:1000000 BP 10000
```

## Library

```rust,no_run
use straw::{straw_stream, MatrixType, Normalization, Unit};

straw_stream(
    MatrixType::Observed,
    Normalization::none(),
    "sample.hic",
    "chr1:0:1000000",
    "chr1:0:1000000",
    Unit::BP,
    10_000,
    |record| println!("{}\t{}\t{}", record.bin_x, record.bin_y, record.counts),
)?;
# Ok::<(), straw::Error>(())
```

Open `HicFile` directly when making repeated queries so the parsed header and
HTTP connection pool are reused. `HicFile::records`, `stream_records`,
`matrix`, `count_records`, and `chromosome_record_counts` correspond to the
C++ library operations.

Header and footer metadata can be inspected without reading contact blocks:

```rust,no_run
use straw::{HicFile, Unit};

let hic = HicFile::open("sample.hic")?;
println!("genome: {}, format: v{}", hic.genome_id(), hic.version());

for chromosome in hic.chromosomes() {
    println!("{}\t{}", chromosome.name, chromosome.length);
}

println!("BP resolutions: {:?}", hic.bp_resolutions());
println!("FRAG resolutions: {:?}", hic.fragment_resolutions());
println!("all normalizations: {:?}", hic.normalizations()?);
println!(
    "chr1 at 10kb: {:?}",
    hic.normalizations_for("chr1", Unit::BP, 10_000)?,
);

for entry in hic.normalization_entries()? {
    println!(
        "{} {} {} {}",
        entry.normalization,
        entry.chromosome.name,
        entry.unit,
        entry.resolution,
    );
}
# Ok::<(), straw::Error>(())
```

A complete metadata example is also available:

```bash
cargo run --release --example metadata -- sample.hic
```

## Compatibility testing

The implementation is tested by comparing CLI output with `C++/build/straw`
on the repository's large `.hic` fixture. Run the standard Rust checks with:

```bash
cargo test
cargo clippy --all-targets -- -D warnings
```
