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
            return SeqReader {stream: Box::new(FastqStream::<GzDecoder<File>>::new(filename))};
        } else if filename.ends_with("fastq"){
            return SeqReader {stream: Box::new(FastqStream::<File>::new(filename))};
        } else if filename.ends_with("fna.gz"){
            return SeqReader {stream: Box::new(FastaStream::<GzDecoder<File>>::new(filename))};
        } else if filename.ends_with("fna"){
            return SeqReader {stream: Box::new(FastaStream::<File>::new(filename))};
        } else{
            panic!("Could not determine the format of file {}", filename);
        }
    }
    
    // New from stdin
    pub fn new_from_stdin(fastq: bool, gzipped: bool) -> SeqReader{
        return SeqReader {stream: Box::new(FastqStream::<std::io::Stdin>::new())};
    /*
        if fastq && gzipped {
            return SeqReader {stream: Box::new(FastqStream::<GzDecoder<std::io::Stdin>>::new(io::stdin()))};
        } else if fastq && !gzipped{
            return SeqReader {stream: Box::new(FastqStream::<std::io::Stdin>::new(io::stdin()))};
        } else if !fastq && gzipped{
            return SeqReader {stream: Box::new(FastaStream::<GzDecoder<std::io::Stdin>>::new(io::stdin()))};
        } else if !fastq && !gzipped{
            return SeqReader {stream: Box::new(FastaStream::<std::io::Stdin>::new(io::stdin()))};
        } else{
            panic!("This line should never be reached");
        }*/
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
        .get_matches();
        
    if let Some(infile) = matches.get_one::<String>("input"){
        println!("inputfile: {:?}", infile);    
    } else{
        println!("No input file");
    }

    println!("fasta flag: {:?}", matches.get_flag("fasta"));

    let args: Vec<String> = env::args().collect();
    
    dbg!(args.len());
    if args.len() == 1{
        // From Stdin
        let mut reader = SeqReader::new_from_stdin(false, true);
        reader.read_all();
    } else{
        // From file
        let filename: String = args[1].clone();
        let mut reader = SeqReader::new(&filename);
        reader.read_all();
    }


}
