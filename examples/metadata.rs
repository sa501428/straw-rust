use straw::{HicFile, Result, Unit};

fn main() -> Result<()> {
    let path = std::env::args()
        .nth(1)
        .ok_or_else(|| straw::Error::Argument("usage: metadata <file.hic>".into()))?;
    let hic = HicFile::open(path)?;

    println!("genome\t{}", hic.genome_id());
    println!("version\t{}", hic.version());
    println!("bp_resolutions\t{:?}", hic.bp_resolutions());
    println!("fragment_resolutions\t{:?}", hic.fragment_resolutions());
    println!("normalizations\t{:?}", hic.normalizations()?);
    if hic.chromosome("chr1").is_some() && hic.bp_resolutions().contains(&10_000) {
        println!(
            "chr1_10kb_normalizations\t{:?}",
            hic.normalizations_for("chr1", Unit::BP, 10_000)?
        );
    }

    for (key, value) in hic.attributes() {
        println!("attribute\t{key}\t{} bytes", value.len());
    }
    for chromosome in hic.chromosomes() {
        println!(
            "chromosome\t{}\t{}\t{}",
            chromosome.index, chromosome.name, chromosome.length
        );
    }
    Ok(())
}
