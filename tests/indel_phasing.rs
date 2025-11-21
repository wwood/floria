use assert_cmd::Command;
use rust_htslib::bam;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Header;
use rust_htslib::bam::Writer;
use rust_htslib::bam::{index, Format};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;
use std::time::SystemTime;

fn make_work_dir() -> PathBuf {
    let suffix = SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("target")
        .join(format!("indel_phase_{}", suffix));
    let _ = fs::remove_dir_all(&base);
    fs::create_dir_all(&base).unwrap();
    base
}

fn write_reference(work_dir: &PathBuf, sequence: &str) -> PathBuf {
    let fasta_path = work_dir.join("ref.fa");
    let mut fasta = File::create(&fasta_path).unwrap();
    writeln!(fasta, ">chr1").unwrap();
    writeln!(fasta, "{}", sequence).unwrap();

    let fai_path = work_dir.join("ref.fa.fai");
    let offset = ">chr1\n".len();
    let line_bases = sequence.len();
    let line_bytes = line_bases + 1;
    let mut fai = File::create(fai_path).unwrap();
    writeln!(
        fai,
        "chr1\t{}\t{}\t{}\t{}",
        sequence.len(),
        offset,
        line_bases,
        line_bytes
    )
    .unwrap();

    fasta_path
}

fn write_vcf(work_dir: &PathBuf) -> PathBuf {
    let vcf_path = work_dir.join("vars.vcf");
    let mut vcf = File::create(&vcf_path).unwrap();
    writeln!(vcf, "##fileformat=VCFv4.2").unwrap();
    writeln!(vcf, "##contig=<ID=chr1,length=40>").unwrap();
    writeln!(
        vcf,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample"
    )
    .unwrap();
    writeln!(vcf, "chr1\t5\t.\tA\tATT\t.\tPASS\t.\tGT\t1/1").unwrap();
    writeln!(vcf, "chr1\t13\t.\tACG\tA\t.\tPASS\t.\tGT\t1/1").unwrap();
    writeln!(vcf, "chr1\t25\t.\tA\tG\t.\tPASS\t.\tGT\t1/1").unwrap();
    vcf_path
}

fn make_header() -> Header {
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", &"1.5");
    hd.push_tag(b"SO", &"coordinate");
    header.push_record(&hd);

    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", &"chr1");
    sq.push_tag(b"LN", &40);
    header.push_record(&sq);
    header
}

fn hap1_cigar() -> bam::record::CigarString {
    vec![
        Cigar::Match(5),
        Cigar::Ins(2),
        Cigar::Match(8),
        Cigar::Del(2),
        Cigar::Match(25),
    ]
    .into()
}

fn hap1_sequence(reference: &[u8]) -> Vec<u8> {
    let mut seq = Vec::new();
    seq.extend_from_slice(&reference[0..5]);
    seq.extend_from_slice(b"TT");
    seq.extend_from_slice(&reference[5..13]);
    seq.extend_from_slice(&reference[15..]);
    seq
}

fn hap2_cigar() -> bam::record::CigarString {
    vec![Cigar::Match(40)].into()
}

fn hap2_sequence(reference: &[u8]) -> Vec<u8> {
    let mut seq = reference.to_vec();
    seq[24] = b'G';
    seq
}

fn write_bam(work_dir: &PathBuf, reference: &[u8]) -> PathBuf {
    let bam_path = work_dir.join("reads.bam");
    let mut writer = Writer::from_path(&bam_path, &make_header(), Format::Bam).unwrap();

    let hap1_seq = hap1_sequence(reference);
    let hap2_seq = hap2_sequence(reference);
    let h1_cigar = hap1_cigar();
    let h2_cigar = hap2_cigar();

    for idx in 0..10 {
        let mut record = bam::Record::new();
        record.set(
            format!("hap1_read_{}", idx).as_bytes(),
            Some(&h1_cigar),
            &hap1_seq,
            &vec![60; hap1_seq.len()],
        );
        record.set_tid(0);
        record.set_pos(0);
        record.set_mapq(60);
        writer.write(&record).unwrap();
    }

    for idx in 0..10 {
        let mut record = bam::Record::new();
        record.set(
            format!("hap2_read_{}", idx).as_bytes(),
            Some(&h2_cigar),
            &hap2_seq,
            &vec![60; hap2_seq.len()],
        );
        record.set_tid(0);
        record.set_pos(0);
        record.set_mapq(60);
        writer.write(&record).unwrap();
    }

    drop(writer);
    index::build(&bam_path, None, index::Type::Bai, 1).unwrap();
    bam_path
}

#[test]
fn floria_phases_reads_with_indels() {
    let work_dir = make_work_dir();
    let reference_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let ref_path = write_reference(&work_dir, reference_seq);
    let vcf_path = write_vcf(&work_dir);
    let bam_path = write_bam(&work_dir, reference_seq.as_bytes());
    let out_dir = work_dir.join("out");

    Command::cargo_bin("floria")
        .unwrap()
        .args([
            "-b",
            bam_path.to_str().unwrap(),
            "-v",
            vcf_path.to_str().unwrap(),
            "-r",
            ref_path.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "--snp-count-filter",
            "1",
        ])
        .assert()
        .success();

    let haplosets_path = out_dir.join("chr1").join("chr1.haplosets");
    let haplosets = fs::read_to_string(&haplosets_path).unwrap();

    let mut current_hap: Option<usize> = None;
    let mut assignments: HashMap<String, usize> = HashMap::new();
    for line in haplosets.lines() {
        if line.starts_with(">HAP") {
            if let Some(rest) = line.strip_prefix(">HAP") {
                let hap_id: usize = rest
                    .split_whitespace()
                    .next()
                    .and_then(|tok| tok.split('.').next())
                    .unwrap()
                    .parse()
                    .unwrap();
                current_hap = Some(hap_id);
            }
        } else if line.starts_with('#') || line.trim().is_empty() {
            continue;
        } else if let Some(hap) = current_hap {
            let read_id = line.split_whitespace().next().unwrap().to_string();
            assignments.insert(read_id, hap);
        }
    }

    let mut hap1_labels = HashSet::new();
    let mut hap2_labels = HashSet::new();
    for idx in 0..10 {
        let hap = assignments.get(&format!("hap1_read_{}", idx)).copied();
        hap1_labels.insert(hap);
        let hap = assignments.get(&format!("hap2_read_{}", idx)).copied();
        hap2_labels.insert(hap);
    }

    assert_eq!(hap1_labels.len(), 1, "haplotype 1 reads clustered together");
    assert_eq!(hap2_labels.len(), 1, "haplotype 2 reads clustered together");
    assert_ne!(hap1_labels, hap2_labels, "haplotypes separated");
}
