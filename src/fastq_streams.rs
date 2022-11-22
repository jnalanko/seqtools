use std::io;
use flate2::read::GzDecoder;
use seq_io::fastq::{Reader,Record};
use std::fs::File;

pub trait SeqStream{
    fn read_all(&mut self);
}

pub struct GzippedFastqStream{
    file_gz_reader: seq_io::fastq::Reader<GzDecoder<std::fs::File>>
}

pub struct RegularFastqStream{
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
    pub fn new(filename: &String) -> GzippedFastqStream{
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
    pub fn new(filename: &String) -> RegularFastqStream{
        return RegularFastqStream{
            file_reader: Reader::new(File::open(&filename).unwrap())
        };
    }
}
