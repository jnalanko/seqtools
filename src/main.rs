extern crate flate2;

use crate::fastq_streams::{SeqStream,FastqStream,FastaStream};

mod fastq_streams;

use std::io;
use flate2::read::GzDecoder;
use std::fs::File;
use std::env;
use std::path::Path;
use std::ffi::OsStr;

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

//    let mut reader = Reader::from_path("/home/niklas/data/SRR19749835_prefix.fastq").unwrap();

/*
    let stdin = io::stdin();
    let d = GzDecoder::new(stdin);
    let mut reader = Reader::new(d);

    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");
        println!("{}", record.id().unwrap());
    }
*/
}
