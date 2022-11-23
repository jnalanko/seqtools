extern crate flate2;

use crate::fastq_streams::{SeqStream,FastqStream,FastaStream};

mod fastq_streams;

use std::io;
use flate2::read::GzDecoder;
use std::fs::File;
use std::env;
use std::path::Path;
use std::ffi::OsStr;
use clap::{Arg,Command,ArgAction};

struct SeqReader {
    stream: Box<dyn SeqStream>
}

impl SeqReader{

    // New from file
    pub fn new(filename: &String) -> SeqReader{
        if filename.ends_with("fastq.gz"){
            return SeqReader {stream: Box::new(FastqStream::<GzDecoder<File>>::new_from_file(filename))};
        } else if filename.ends_with("fastq"){
            return SeqReader {stream: Box::new(FastqStream::<File>::new_from_file(filename))};
        } else if filename.ends_with("fna.gz"){
            return SeqReader {stream: Box::new(FastaStream::<GzDecoder<File>>::new_from_file(filename))};
        } else if filename.ends_with("fna"){
            return SeqReader {stream: Box::new(FastaStream::<File>::new_from_file(filename))};
        } else{
            panic!("Could not determine the format of file {}", filename);
        }
    }
    
    // New from stdin
    pub fn new_from_stdin(fastq: bool, gzipped: bool) -> SeqReader{
        if fastq && gzipped {
            return SeqReader {stream: Box::new(FastqStream::<GzDecoder<std::io::Stdin>>::new(GzDecoder::new(io::stdin())))};
        } else if fastq && !gzipped{
            return SeqReader {stream: Box::new(FastqStream::<std::io::Stdin>::new(io::stdin()))};
        } else if !fastq && gzipped{
            return SeqReader {stream: Box::new(FastaStream::<GzDecoder<std::io::Stdin>>::new(GzDecoder::new(io::stdin())))};
        } else if !fastq && !gzipped{
            return SeqReader {stream: Box::new(FastaStream::<std::io::Stdin>::new(io::stdin()))};
        } else{
            panic!("This line should never be reached");
        }
    }
    
    pub fn read_all(&mut self){
        self.stream.read_all();
    }
}


fn main(){
    let matches = Command::new("Fasta/fastq parsing")
        .version("0.1.0")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .about("Fasta/fastq parsing")
        .arg(Arg::new("input")
                 .short('i')
                 .long("input")
                 .help("Input filename"))
        .arg(Arg::new("fasta")
                 .short('a')
                 .long("fasta")
                 .action(ArgAction::SetTrue)
                 .help("Parse in fasta format"))
        .arg(Arg::new("fastq")
                 .short('q')
                 .long("fastq")
                 .action(ArgAction::SetTrue)
                 .help("Parse in fastq format"))
        .arg(Arg::new("gzip")
                 .short('g')
                 .long("gzip")
                 .action(ArgAction::SetTrue)
                 .help("Read gzipped input"))
        .get_matches();
        
    if let Some(infile) = matches.get_one::<String>("input"){
        let mut reader = SeqReader::new(&infile);
        reader.read_all();
    } else{
        let is_fasta = matches.get_flag("fasta");
        let is_fastq = matches.get_flag("fastq");
        let is_gzip = matches.get_flag("gzip");
        if is_fasta && is_fastq {
            println!("Error: can't give both fasta and fastq flags.");
            std::process::exit(-1)
        }
        if !is_fasta && !is_fastq {
            println!("Error: must give --fasta or --fastq and possibly --gzip if reading from stdin.");
            std::process::exit(-1)
        }
        let mut reader = SeqReader::new_from_stdin(is_fastq, is_gzip);
        reader.read_all();
    }
}
