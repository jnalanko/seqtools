extern crate flate2;

use std::io;
use flate2::read::GzDecoder;
use seq_io::fastq::{Reader,Record};

struct SeqReader {
    gzipped: bool,
    reader: seq_io::fastq::Reader<std::io::Stdin>,
    gz_reader: seq_io::fastq::Reader<GzDecoder<std::io::Stdin>>
}

impl SeqReader{

    pub fn new(gzipped: bool) -> SeqReader{
        SeqReader {gzipped: gzipped, reader: Reader::new(io::stdin()), gz_reader: Reader::new(GzDecoder::new(io::stdin()))}
    }
    
    fn read_all(&mut self){
        if self.gzipped {
            while let Some(record) = self.reader.next() {
                let record = record.expect("Error reading record");
                println!("{}", record.id().unwrap());
            }
        } else{
            while let Some(record) = self.gz_reader.next() {
                let record = record.expect("Error reading record");
                println!("{}", record.id().unwrap());
            }
        }
    }

}


fn main(){

    let mut reader = SeqReader::new(true);
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
