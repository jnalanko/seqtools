extern crate flate2;

use std::io;
use flate2::read::GzDecoder;
use seq_io::fastq::{Reader,Record};
use std::fs::File;
use std::env;

struct SeqReader {
    gzipped: bool,
    from_file: bool,
    stdin_reader: seq_io::fastq::Reader<std::io::Stdin>,
    stdin_gz_reader: seq_io::fastq::Reader<GzDecoder<std::io::Stdin>>,
    file_reader: seq_io::fastq::Reader<std::fs::File>,
    file_gz_reader: seq_io::fastq::Reader<GzDecoder<std::fs::File>>
}

impl SeqReader{

    pub fn new(filename: &String, gzipped: bool) -> SeqReader{
        println!("{} {}", gzipped, filename.chars().count() > 0);
        SeqReader {gzipped: gzipped, 
                   from_file: (filename.chars().count() > 0),
                   stdin_reader: Reader::new(io::stdin()),
                   stdin_gz_reader: Reader::new(GzDecoder::new(io::stdin())),
                   file_reader: Reader::new(File::open(&filename).unwrap()),
                   file_gz_reader: Reader::new(GzDecoder::new(File::open(&filename).unwrap()))
        }
    }
    
    fn read_all(&mut self){
        
        if self.gzipped && self.from_file {
            println!("A");
            while let Some(record) = self.file_gz_reader.next() {
                let record = record.expect("Error reading record");
                println!("{}", record.id().unwrap());
            }
        }
        
        
        if self.gzipped && !self.from_file {
            println!("B");
            while let Some(record) = self.stdin_gz_reader.next() {
                let record = record.expect("Error reading record");
                println!("{}", record.id().unwrap());
            }
        }
        
        if !self.gzipped && self.from_file {
            println!("C");
            while let Some(record) = self.file_reader.next() {
                let record = record.expect("Error reading record");
                println!("{}", record.id().unwrap());
            }
        }
        
        if !self.gzipped && !self.from_file {
            println!("D");
            while let Some(record) = self.stdin_reader.next() {
                let record = record.expect("Error reading record");
                println!("{}", record.id().unwrap());
            }
        }
    }
}


fn main(){

    let args: Vec<String> = env::args().collect();
    let filename: String = args[1].clone();
    let n = filename.chars().count();
    println!("{}", n);
    let gzipped = &filename[n-3..n] == ".gz";
    println!("{}", filename);
    println!("{}", gzipped);
    
    let mut file_test = match File::open(&filename){
        Err(why) => panic!("couldn't read: {}", why),
        Ok(file_test) => file_test
    };
    
    let mut gz_file_test = GzDecoder::new(file_test);
    
    println!("Before");
    let mut reader = SeqReader::new(&filename, gzipped);
    println!("After");    
    
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
