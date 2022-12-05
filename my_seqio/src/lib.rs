
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
pub mod record;
mod tests;

use crate::reader::FastXReader;

#[derive(Copy, Clone)]
pub enum FileType{
    FASTA,
    FASTQ,
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
