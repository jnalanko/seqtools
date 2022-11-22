extern crate flate2;

use std::io;
use flate2::read::GzDecoder;
use seq_io::fastq::{Reader,Record};
use std::fs::File;
use std::env;

trait SeqStream{

    fn read_all(&mut self);

}

struct GzippedFastqStream{
    file_gz_reader: seq_io::fastq::Reader<GzDecoder<std::fs::File>>
}

struct RegularFastqStream{
    file_reader: seq_io::fastq::Reader<std::fs::File>
}

impl SeqStream for GzippedFastqStream{

    fn read_all(&mut self){
        while let Some(record) = self.file_gz_reader.next() {
            let record = record.expect("Error reading record");
            println!("{}", record.id().unwrap());
        }
    }
}

impl GzippedFastqStream{
    fn new(filename: &String) -> GzippedFastqStream{
        return GzippedFastqStream{
            file_gz_reader: Reader::new(GzDecoder::new(File::open(&filename).unwrap()))
        };
    }
}

impl SeqStream for RegularFastqStream{

    fn read_all(&mut self){
        while let Some(record) = self.file_reader.next() {
            let record = record.expect("Error reading record");
            println!("{}", record.id().unwrap());
        }
    }
}

impl RegularFastqStream{
    fn new(filename: &String) -> RegularFastqStream{
        return RegularFastqStream{
            file_reader: Reader::new(File::open(&filename).unwrap())
        };
    }
}

struct SeqReader {
    stream: Box<dyn SeqStream>
}

impl SeqReader{

    pub fn new(filename: &String) -> SeqReader{
        let n = filename.chars().count();
        let gzipped = &filename[n-3..n] == ".gz";
        if gzipped{
            return SeqReader {stream: Box::new(GzippedFastqStream::new(filename))};
        } else {
            return SeqReader {stream: Box::new(RegularFastqStream::new(filename))};
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
