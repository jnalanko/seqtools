
use std::io;
use std::io::BufReader;
use std::io::BufRead;
use std::io::BufWriter;
use std::io::Write;
use std::fmt;
use std::str;
use std::fs::File;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;

pub mod reader;
pub mod writer;

use crate::reader::FastXReader;

#[derive(Copy, Clone)]
pub enum FileType{
    FASTA,
    FASTQ,
}

pub trait Record{
    fn head(&self) -> &[u8];
    fn seq(&self) -> &[u8];
    fn qual(&self) -> Option<&[u8]>;
}

#[derive(Debug)]
pub struct SeqRecord<'a>{
    pub head: &'a [u8],    
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>, // If FASTA, this is None
}

#[derive(Debug)]
pub struct OwnedSeqRecord{
    pub head: Vec<u8>,    
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>, // If FASTA, this is None
}

impl<'a> Record for SeqRecord<'a>{
    fn head(&self) -> &[u8]{self.head}
    fn seq(&self) -> &[u8]{self.seq}
    fn qual(&self) -> Option<&[u8]>{self.qual}
}

impl<'a> Record for OwnedSeqRecord{
    fn head(&self) -> &[u8]{self.head.as_slice()}
    fn seq(&self) -> &[u8]{self.seq.as_slice()}
    fn qual(&self) -> Option<&[u8]>{
        match &self.qual{
            Some(q) => return Some(q.as_slice()),
            None => None,
        }
    }
}

impl<'a> SeqRecord<'a>{
    pub fn to_owned(&self) -> OwnedSeqRecord{
        OwnedSeqRecord { 
            head: self.head.to_vec(), 
            seq: self.seq.to_vec(), 
            qual: match self.qual {
                Some(q) => Some(q.to_vec()), 
                None => None
            }
        }
    }
}

impl<'a> fmt::Display for SeqRecord<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
               "SeqRecord{{ \n  Head: {}\n  Seq:  {}\n  Qual: {}\n}}", 
               str::from_utf8(self.head).unwrap(),
               str::from_utf8(self.seq).unwrap(),
               match self.qual{
                   Some(q) => str::from_utf8(q).unwrap(),
                   None => "", // No quality values
               }
               
        )
    }
}

// Returns (file type, is_gzipped)
fn figure_out_file_format(filename: &str) -> (FileType, bool){
    let is_gzipped = filename.ends_with(".gz");
    let filename = if is_gzipped{
        &filename[0 .. filename.len()-3] // Drop the .gz suffix
    } else {&filename};
    let fasta_extensions = vec![".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa"];
    let fastq_extensions = vec![".fastq", ".fq"];
    if fasta_extensions.iter().any(|&suffix| filename.ends_with(suffix)){
        return (FileType::FASTA, is_gzipped);
    } else if fastq_extensions.iter().any(|&suffix| filename.ends_with(suffix)){
        return (FileType::FASTQ, is_gzipped);
    } else{
        panic!("Unkown file extension: {}", filename);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{cmp::min, process::Output};

    #[test]
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

    #[test]
    fn test_figure_out_file_format(){
        assert!(match figure_out_file_format("aa.fna") {(FileType::FASTA,false) => true, _ => false});
        assert!(match figure_out_file_format("aa.fq") {(FileType::FASTQ,false) => true, _ => false});
        assert!(match figure_out_file_format("bbb.fna.gz") {(FileType::FASTA,true) => true, _ => false});
        assert!(match figure_out_file_format("cc.fna.gz") {(FileType::FASTA,true) => true, _ => false});
        assert!(match figure_out_file_format(".fna.gz") {(FileType::FASTA,true) => true, _ => false});
        assert!(match figure_out_file_format(".fasta") {(FileType::FASTA,false) => true, _ => false});
        assert!(match figure_out_file_format(".fq") {(FileType::FASTQ,false) => true, _ => false});
    }
}
