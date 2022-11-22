extern crate flate2;

use crate::fastq_streams::{SeqStream,GzippedFastqStream,RegularFastqStream,TemplateTest};

mod fastq_streams;

use std::io;
use flate2::read::GzDecoder;
use seq_io::fastq::{Reader,Record};
use std::fs::File;
use std::env;

struct SeqReader {
    stream: Box<dyn SeqStream>
}

impl SeqReader{

    pub fn new(filename: &String) -> SeqReader{
        let n = filename.chars().count();
        let gzipped = &filename[n-3..n] == ".gz";
        if gzipped{
            return SeqReader {stream: Box::new(TemplateTest::<GzDecoder<File>>::new(filename))};
        } else {
            return SeqReader {stream: Box::new(TemplateTest::<File>::new(filename))};
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
