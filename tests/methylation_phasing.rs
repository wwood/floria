use assert_cmd::Command;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Aux, Cigar, CigarString, Record};
use rust_htslib::bam::{self, Format, Header, Writer};
use std::fs;
use std::io::Write;
use std::path::Path;
use tempfile::tempdir;

fn write_reference(dir: &Path, seq: &str) -> (String, String) {
    let fasta_path = dir.join("ref.fa");
    let mut fasta = fs::File::create(&fasta_path).unwrap();
    writeln!(fasta, ">ref\n{}", seq).unwrap();
    let fasta_str = fasta_path.to_string_lossy().to_string();
    let fai_path = fasta_path.with_extension("fa.fai");
    let mut fai = fs::File::create(&fai_path).unwrap();
    // For this test data we write a single-line fasta; offsets are straightforward.
    let offset = 5; // length of ">ref\n"
    let line_bases = seq.len();
    let line_width = line_bases + 1; // include trailing newline
    writeln!(
        fai,
        "ref\t{}\t{}\t{}\t{}",
        line_bases, offset, line_bases, line_width
    )
    .unwrap();
    (fasta_str.clone(), fai_path.to_string_lossy().to_string())
}

fn build_record(name: &str, seq: &[u8], pos: i64, mm_tag: Option<&str>) -> Record {
    let mut record = Record::new();
    record.set_qname(name.as_bytes());
    record.set_tid(0);
    record.set_pos(pos);
    record.set_mapq(60);
    record.set_flags(0);
    let cigar = CigarString::from(vec![Cigar::Match(seq.len() as u32)]);
    record.set(name.as_bytes(), Some(&cigar), seq, &vec![30u8; seq.len()]);
    if let Some(tag) = mm_tag {
        record
            .push_aux(b"MM", Aux::String(tag))
            .expect("failed to add MM tag");
    }
    record
}

fn write_bam(dir: &Path, reads: Vec<Record>) -> String {
    let bam_path = dir.join("reads.bam");
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", &"1.6");
    hd.push_tag(b"SO", &"coordinate");
    header.push_record(&hd);
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", &"ref");
    sq.push_tag(b"LN", &40);
    header.push_record(&sq);
    let mut writer = Writer::from_path(&bam_path, &header, Format::Bam).unwrap();
    for mut rec in reads {
        writer.write(&mut rec).unwrap();
    }
    drop(writer);
    bam::index::build(bam_path.to_str().unwrap(), None, bam::index::Type::Bai, 1).unwrap();
    bam_path.to_string_lossy().to_string()
}

fn write_vcf(dir: &Path, records: &[(&str, u32, char, char)]) -> String {
    let vcf_path = dir.join("variants.vcf");
    let mut vcf = fs::File::create(&vcf_path).unwrap();
    writeln!(vcf, "##fileformat=VCFv4.2").unwrap();
    writeln!(vcf, "##contig=<ID=ref,length=40>").unwrap();
    writeln!(
        vcf,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    )
    .unwrap();
    for (chrom, pos, ref_base, alt_base) in records {
        writeln!(
            vcf,
            "{}\t{}\t.\t{}\t{}\t50\tPASS\t.\tGT\t0/1",
            chrom, pos, ref_base, alt_base
        )
        .unwrap();
    }
    vcf_path.to_string_lossy().to_string()
}

fn write_bedmethyl(dir: &Path, records: &[(&str, u32, &str)]) -> String {
    let bed_path = dir.join("sites.bed");
    let mut bed = fs::File::create(&bed_path).unwrap();
    for (chrom, start, code) in records {
        writeln!(bed, "{}\t{}\t{}\t{}", chrom, start, start + 1, code).unwrap();
    }
    bed_path.to_string_lossy().to_string()
}

fn run_floria(reference: &str, vcf: &str, bed: &str, bam: &str, out_dir: &Path) {
    Command::cargo_bin("floria")
        .unwrap()
        .args([
            "-r",
            reference,
            "-v",
            vcf,
            "--bedmethyl",
            bed,
            "-b",
            bam,
            "-o",
            out_dir.to_str().unwrap(),
            "--snp-count-filter",
            "0",
            "--block-length",
            "10",
            "--snp-density",
            "0.0",
            "-t",
            "1",
            "--max-ploidy",
            "2",
            "--no-stop-heuristic",
        ])
        .assert()
        .success();
}

fn read_haplosets(out_dir: &Path) -> Vec<Vec<String>> {
    let mut haplosets = vec![];
    for entry in walkdir::WalkDir::new(out_dir) {
        let entry = entry.unwrap();
        if entry.file_name() == "ref.haplosets" {
            let content = fs::read_to_string(entry.path()).unwrap();
            let mut current = vec![];
            for line in content.lines() {
                if line.starts_with('>') {
                    if !current.is_empty() {
                        haplosets.push(current);
                        current = vec![];
                    }
                } else if !line.is_empty() && !line.starts_with('#') {
                    let cols: Vec<&str> = line.split('\t').collect();
                    current.push(cols[0].to_string());
                }
            }
            if !current.is_empty() {
                haplosets.push(current);
            }
        }
    }
    haplosets
}

fn build_sequences(
    template: &[u8],
    methyl_c_positions: &[usize],
    snp_override: Option<(usize, u8)>,
) -> Vec<u8> {
    let mut seq = template.to_vec();
    for pos in methyl_c_positions {
        seq[*pos] = b'C';
    }
    if let Some((pos, base)) = snp_override {
        seq[pos] = base;
    }
    seq
}

#[test]
fn phases_snp_and_methylated_reads() {
    let dir = tempdir().unwrap();
    let reference_seq = {
        let mut seq = vec![b'A'; 40];
        seq[9] = b'C';
        seq[19] = b'C';
        String::from_utf8(seq.clone()).unwrap()
    };
    let (fasta_path, _) = write_reference(dir.path(), &reference_seq);

    let template = reference_seq.as_bytes().to_vec();
    let hap1_seq = build_sequences(&template, &[9, 19], None);
    let hap2_seq = build_sequences(&template, &[], Some((29, b'G')));

    let mut reads = vec![];
    for i in 0..10 {
        reads.push(build_record(
            &format!("hap1_{}", i),
            &hap1_seq,
            0,
            Some("C+m,0,0;"),
        ));
        reads.push(build_record(&format!("hap2_{}", i), &hap2_seq, 0, None));
    }
    let bam_path = write_bam(dir.path(), reads);

    let vcf_path = write_vcf(dir.path(), &[("ref", 30, 'A', 'G')]);
    let bed_path = write_bedmethyl(dir.path(), &[("ref", 9, "m"), ("ref", 19, "m")]);

    let out_dir = dir.path().join("out1");
    run_floria(&fasta_path, &vcf_path, &bed_path, &bam_path, &out_dir);

    let haplosets = read_haplosets(&out_dir);
    assert_eq!(haplosets.len(), 2);
    let mut hap1_count = 0;
    let mut hap2_count = 0;
    for group in haplosets {
        if group.iter().all(|id| id.starts_with("hap1_")) {
            hap1_count = group.len();
        } else if group.iter().all(|id| id.starts_with("hap2_")) {
            hap2_count = group.len();
        }
    }
    assert_eq!(hap1_count, 10);
    assert_eq!(hap2_count, 10);
}

#[test]
fn phases_methylation_per_haplotype() {
    let dir = tempdir().unwrap();
    let reference_seq = {
        let mut seq = vec![b'A'; 40];
        seq[9] = b'C';
        seq[19] = b'C';
        seq[29] = b'A';
        seq[34] = b'C';
        String::from_utf8(seq.clone()).unwrap()
    };
    let (fasta_path, _) = write_reference(dir.path(), &reference_seq);

    let template = reference_seq.as_bytes().to_vec();
    let hap1_seq = build_sequences(&template, &[9], Some((29, b'G')));
    let hap2_seq = build_sequences(&template, &[], Some((34, b'T')));

    let mut reads = vec![];
    for i in 0..10 {
        reads.push(build_record(
            &format!("hap1b_{}", i),
            &hap1_seq,
            0,
            Some("C+21839,0;"),
        ));
        reads.push(build_record(
            &format!("hap2b_{}", i),
            &hap2_seq,
            0,
            Some("A+a,27;"),
        ));
    }
    let bam_path = write_bam(dir.path(), reads);

    let vcf_path = write_vcf(dir.path(), &[("ref", 30, 'A', 'G'), ("ref", 35, 'A', 'T')]);
    let bed_path = write_bedmethyl(dir.path(), &[("ref", 9, "21839"), ("ref", 29, "a")]);

    let out_dir = dir.path().join("out2");
    run_floria(&fasta_path, &vcf_path, &bed_path, &bam_path, &out_dir);

    let haplosets = read_haplosets(&out_dir);
    assert_eq!(haplosets.len(), 2);
    let mut hap1_count = 0;
    let mut hap2_count = 0;
    for group in haplosets {
        if group.iter().all(|id| id.starts_with("hap1b_")) {
            hap1_count = group.len();
        } else if group.iter().all(|id| id.starts_with("hap2b_")) {
            hap2_count = group.len();
        }
    }
    assert_eq!(hap1_count, 10);
    assert_eq!(hap2_count, 10);
}
