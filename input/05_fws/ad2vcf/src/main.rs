use npyz::NpyFile;
use npyz::Order;
use rust_htslib::bcf::header::Header;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{Format, Writer};
use std::io::BufReader;

fn to_npy_reader(filename: &str) -> NpyFile<BufReader<std::fs::File>> {
    std::fs::File::open(filename)
        .map(BufReader::new)
        .map(NpyFile::new)
        .unwrap()
        .unwrap()
}

fn get_positions(filename: &str) -> Vec<i32> {
    to_npy_reader(filename).into_vec::<i32>().unwrap()
}

/// numpy array that contains object type should be converted to '|SN' type before save to a npy
/// file where `N` is the maxinum length of the elements including the null-terminated byte
fn get_alleles(filename: &str) -> Vec<String> {
    to_npy_reader(filename).into_vec::<String>().unwrap()
}
/// numpy array that contains object type should be converted to '|S' type before save to a npy
fn get_chroms(filename: &str) -> Vec<String> {
    to_npy_reader(filename).into_vec::<String>().unwrap()
}
fn get_samples(filename: &str) -> Vec<String> {
    to_npy_reader(filename).into_vec::<String>().unwrap()
}

fn main() {
    let chrom = get_chroms("./chr.npy");
    let pos = get_positions("./pos.npy");
    let alleles = get_alleles("alleles.npy");
    let samples = get_samples("./samples.npy");
    assert_eq!(pos.len(), chrom.len());
    assert_eq!(
        pos.len() * 2,
        alleles.len(),
        "site should be biallelic only"
    );

    let ad = to_npy_reader("../pf7_ad.npy");
    let gt = to_npy_reader("../pf7_gt.npy");

    let shape = ad.shape().to_vec();
    let order = ad.order();
    assert_eq!(order, Order::C);
    assert_eq!(shape[0] as usize, pos.len());
    assert_eq!(shape[1] as usize, samples.len());
    assert_eq!(shape[2] as usize, 2, "site should be biallelic only");

    let shape = gt.shape().to_vec();
    let order = gt.order();
    assert_eq!(order, Order::C);
    assert_eq!(shape[0] as usize, pos.len());
    assert_eq!(shape[1] as usize, samples.len());
    assert_eq!(shape[2] as usize, 2, "site should be biallelic only");

    let nsam = samples.len();
    let nsite = pos.len();

    // iter sites
    let mut ad_iter = ad.data::<i16>().unwrap().into_iter();
    let mut pos_iter = pos.into_iter();
    let mut chr_iter = chrom.into_iter();
    let mut alleles_iter = alleles.into_iter();
    let mut gt_iter = gt.data::<i8>().unwrap().into_iter();

    // vcf header
    let lines = vec![
        r#"##fileformat=VCFv4.2"#,
        r#"##FILTER=<ID=PASS,Description="All filters passed">"#,
        r#"##filedate=20230829"#,
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">"#,
        r#"##contig=<ID=Pf3D7_01_v3,length=640851,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_02_v3,length=947102,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_03_v3,length=1067971,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_04_v3,length=1200490,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_05_v3,length=1343557,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_06_v3,length=1418242,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_07_v3,length=1445207,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_08_v3,length=1472805,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_09_v3,length=1541735,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_10_v3,length=1687656,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_11_v3,length=2038340,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_12_v3,length=2271494,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_13_v3,length=2925236,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_14_v3,length=3291936,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf3D7_API_v3,length=34250,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
        r#"##contig=<ID=Pf_M76611,length=5967,assembly=PlasmoDB-44_Pfalciparum3D7_Genome.fasta>"#,
    ];

    let mut header = Header::new();
    for line in lines {
        header.push_record(line.as_bytes());
    }
    for s in samples {
        header.push_sample(s.as_bytes());
    }

    let mut vcf = Writer::from_path("out.bcf", &header, false, Format::Bcf).unwrap();

    let mut record = vcf.empty_record();
    let mut gta_v = Vec::<_>::with_capacity(2 * nsam);
    let mut ad_v = Vec::<_>::with_capacity(2 * nsam);
    let mut n_common_siate = 0;
    for isite in 0..nsite {
        let chrom = chr_iter.next().unwrap();
        let pos = pos_iter.next().unwrap();
        let ref_allele = alleles_iter.next().unwrap();
        let alt_allele = alleles_iter.next().unwrap();

        let mut counter0 = 0;
        let mut counter1 = 0;
        gta_v.clear();
        for _ in 0..(2 * nsam) {
            let gt0: i8 = gt_iter.next().unwrap().unwrap();
            match gt0 {
                0 => {
                    counter0 += 1;
                    gta_v.push(GenotypeAllele::Unphased(0));
                }
                1 => {
                    counter1 += 1;
                    gta_v.push(GenotypeAllele::Unphased(1));
                }
                -1 => {
                    gta_v.push(GenotypeAllele::UnphasedMissing);
                }
                _ => {
                    panic!("unacceptable genotype!")
                }
            }
        }
        let maf = if counter0 < counter1 {
            counter0 as f64 / (counter0 + counter1) as f64
        } else {
            counter1 as f64 / (counter0 + counter1) as f64
        };
        // skip rare variants
        if maf < 0.01 {
            // before end this loop, we need consume some elements from ad
            for _ in 0..(2 * nsam) {
                ad_iter.next();
            }
            continue;
        } else {
            ad_v.clear();
            for _ in 0..(2 * nsam) {
                let ad0 = ad_iter.next().unwrap().unwrap();
                ad_v.push(ad0 as i32);
            }
        }
        record.clear();

        // chrom
        let rid = vcf.header().name2rid(chrom.as_bytes()).unwrap();
        record.set_rid(Some(rid));
        // pos
        record.set_pos(pos as i64 - 1); // 1-based into 0-based
                                        // ref/alt
        record
            .set_alleles(&[ref_allele.as_bytes(), alt_allele.as_bytes()])
            .unwrap();
        // GT
        record.push_genotypes(&gta_v).unwrap();
        // AD
        record.push_format_integer(b"AD", &ad_v[..]).unwrap();

        vcf.write(&record).unwrap();
        n_common_siate += 1;

        println!(
            "isite: {isite}: perc: {:.3}, common variants: {}",
            isite as f32 / nsite as f32,
            n_common_siate
        );
    }
}
