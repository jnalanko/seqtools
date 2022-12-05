use my_seqio::reader::FastXReader;
use my_seqio::writer::FastXWriter;
use my_seqio::record::OwnedSeqRecord;
use my_seqio::FileType;
use my_seqio::figure_out_file_format;
use std::io::BufReader;

use std::{cmp::min, process::Output};

fn fastq() {
    let headers = vec!(
        "SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1",
        "SRR403017.2 HWUSI-EAS108E_0007:3:1:10327:976/1",
        "SRR403017.3 HWUSI-EAS108E_0007:3:1:13569:972/1");
    let seqs = vec!(
        "TTGGACCGGCGCAAGACGGACCAGNGCGAAAGCATTTGCCAAGAANNNN",
        "CAACTTTCTATCTGGCATTCCCTGNGGAGGAAATAGAATGCGCGCNNNN",
        "GATCGGAAGAGCACACGTCTGAACNCCAGTCACTTAGGCATCTCGNNNN",
    );
    let quals = vec!(
        "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQ",
        "RSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~####",
        "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    );

    let n_seqs = headers.len();
    let mut fastq_data: String = "".to_owned();
    for i in 0..n_seqs{
        fastq_data.push_str(format!("@{}\n", headers[i]).as_str());
        fastq_data.push_str(format!("{}\n", seqs[i]).as_str());
        fastq_data.push_str("+\n");
        fastq_data.push_str(format!("{}\n", quals[i]).as_str());
    }

    let input = BufReader::new(fastq_data.as_bytes());
    let mut reader = FastXReader::new(input, FileType::FASTQ);

    let mut owned_records: Vec<OwnedSeqRecord> = vec![];
    let mut seqs_read = 0;
    loop{
        if let Some(record) = reader.next(){
            assert_eq!(record.head, headers[seqs_read].as_bytes());
            assert_eq!(record.seq, seqs[seqs_read].as_bytes());
            assert_eq!(record.qual.unwrap(), quals[seqs_read].as_bytes());
            owned_records.push(record.to_owned());
            seqs_read += 1;
        } else { break };
    }
    assert_eq!(seqs_read, n_seqs);

    // Test writer
    let out_buf: Vec<u8> = vec![];
    let mut writer = FastXWriter::<Vec<u8>>::new(out_buf, FileType::FASTQ);

    for rec in owned_records.iter() {
        writer.write(rec);
    }

    writer.flush();
    let written_data = writer.output.into_inner().unwrap();

    // Read the records back from written data and compare to originals.

    let mut reader2 = FastXReader::new(written_data.as_slice(), FileType::FASTQ);
    let mut seqs_read2 = 0;
    loop{
        if let Some(record) = reader2.next(){
            dbg!(&record);
            assert_eq!(record.head, headers[seqs_read2].as_bytes());
            assert_eq!(record.seq, seqs[seqs_read2].as_bytes());
            assert_eq!(record.qual.unwrap(), quals[seqs_read2].as_bytes());
            seqs_read2 += 1;
        } else { break };
    }
    assert_eq!(seqs_read2, n_seqs);
}

#[test]
fn fasta() {
    let headers: Vec<String> = vec!(
        "SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1".to_owned(),
        "SRR403017.2 HWUSI-EAS108E_0007:3:1:10327:976/1".to_owned(),
        "SRR403017.3 HWUSI-EAS108E_0007:3:1:13569:972/1".to_owned());
    let seqs: Vec<String> = vec!(
        "TTGGACCGGCGCAAGACGGACCAGNGCGAAAGCATTTGCCAAGAANNNN".to_owned(),
        "CAACTTTCTATCTGGCATTCCCTGNGGAGGAAATAGAATGCGCGCNNNN".to_owned(),
        "GATCGGAAGAGCACACGTCTGAACNCCAGTCACTTAGGCATCTCGNNNN".to_owned(),
    );

    fn split_seq_to_lines(seq: &String, line_length: usize) -> Vec<String>{
        let mut i: usize = 0;
        let mut lines = Vec::<String>::new();
        while line_length*i < seq.len(){
            lines.push(seq[line_length*i .. min(line_length*(i+1), seq.len())].to_owned());
            i += 1;
        }
        lines
    }

    let n_seqs = headers.len();
    let mut fasta_data: String = "".to_owned();
    for i in 0..n_seqs{
        fasta_data.push_str(format!(">{}\n", headers[i].as_str()).as_str());
        for line in split_seq_to_lines(&seqs[i], 11){
            // Line length 11 should make it so that the last line has a different
            // length than the other lines.
            fasta_data.push_str(format!("{}\n", line.as_str()).as_str());
        }
    }

    dbg!(&fasta_data);

    let input = BufReader::new(fasta_data.as_bytes());
    let mut reader = FastXReader::new(input, FileType::FASTA);

    let mut owned_records: Vec<OwnedSeqRecord> = vec![];
    let mut seqs_read = 0;
    loop{
        if let Some(record) = reader.next(){
            dbg!(&record);
            assert_eq!(record.head, headers[seqs_read].as_bytes());
            assert_eq!(record.seq, seqs[seqs_read].as_bytes());
            assert_eq!(record.qual, None);
            owned_records.push(record.to_owned());
            seqs_read += 1;
        } else { break };
    }
    assert_eq!(seqs_read, n_seqs);

    // Test writer
    let out_buf: Vec<u8> = vec![];
    let mut writer = FastXWriter::<Vec<u8>>::new(out_buf, FileType::FASTA);

    for rec in owned_records.iter() {
        writer.write(rec);
    }

    writer.flush();
    let written_data = writer.output.into_inner().unwrap();

    // This written data may not exactly equal the original data,
    // because the length of FASTA sequence lines is not fixed.
    // Read the records back from written data and compare to originals.

    let mut reader2 = FastXReader::new(written_data.as_slice(), FileType::FASTA);
    let mut seqs_read2 = 0;
    loop{
        if let Some(record) = reader2.next(){
            dbg!(&record);
            assert_eq!(record.head, headers[seqs_read2].as_bytes());
            assert_eq!(record.seq, seqs[seqs_read2].as_bytes());
            assert_eq!(record.qual, None);
            seqs_read2 += 1;
        } else { break };
    }
    assert_eq!(seqs_read2, n_seqs);

}

fn test_figure_out_file_format(){
    assert!(match figure_out_file_format("aa.fna") {(FileType::FASTA,false) => true, _ => false});
    assert!(match figure_out_file_format("aa.fq") {(FileType::FASTQ,false) => true, _ => false});
    assert!(match figure_out_file_format("bbb.fna.gz") {(FileType::FASTA,true) => true, _ => false});
    assert!(match figure_out_file_format("cc.fna.gz") {(FileType::FASTA,true) => true, _ => false});
    assert!(match figure_out_file_format(".fna.gz") {(FileType::FASTA,true) => true, _ => false});
    assert!(match figure_out_file_format(".fasta") {(FileType::FASTA,false) => true, _ => false});
    assert!(match figure_out_file_format(".fq") {(FileType::FASTQ,false) => true, _ => false});
}
