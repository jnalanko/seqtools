extern crate flate2;

use crate::fastq_streams::{SeqStream,FastqStream};

mod fastq_streams;

use std::io;
use flate2::read::GzDecoder;
use std::fs::File;
use std::env;
use std::path::Path;
use std::ffi::OsStr;

fn get_extension_from_filename(filename: &str) -> Option<&str> {
    Path::new(filename)
        .extension()
        .and_then(OsStr::to_str)
}

struct SeqReader {
    stream: Box<dyn SeqStream>
}

impl SeqReader{

    pub fn new(filename: &String) -> SeqReader{
        let extension = get_extension_from_filename(filename).unwrap();
        let gzipped = extension == "gz";
        println!("{} {}", extension, gzipped);
        if gzipped{
            return SeqReader {stream: Box::new(FastqStream::<GzDecoder<File>>::new(filename))};
        } else {
            return SeqReader {stream: Box::new(FastqStream::<File>::new(filename))};
        }    
    }
    
    pub fn read_all(&mut self){
        self.stream.read_all();
    }
}


fn main(){

    let args: Vec<String> = env::args().collect();
    let filename: String = args[1].clone();
    let n = filename.chars().count();
    let gzipped = &filename[n-3..n] == ".gz";
    
    let mut reader = SeqReader::new(&filename);
    reader.read_all();

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
