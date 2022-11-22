use std::io;
use flate2::read::GzDecoder;
use seq_io::fastq::{Reader,Record};
use std::fs::File;

pub trait SeqStream{
    fn read_all(&mut self);
}

pub struct TemplateTest<T: std::io::Read>{
    reader: seq_io::fastq::Reader<T>
}

impl<T: std::io::Read> SeqStream for TemplateTest<T>{

    fn read_all(&mut self){
        let mut sum: i64 = 0;
        while let Some(record) = self.reader.next() {
            let record = record.expect("Error reading record");
            sum += record.seq().len() as i64;
        }
        println!("{}", sum);
    }
}

impl TemplateTest<File>{
    pub fn new(filename: &String) -> TemplateTest<File>{
        return TemplateTest::<File>{
            reader: Reader::new(File::open(&filename).unwrap())
        };
    }
}

impl TemplateTest<GzDecoder<File>>{
    pub fn new(filename: &String) -> TemplateTest<GzDecoder<File>>{
        return TemplateTest::<GzDecoder<File>>{
            reader: Reader::new(GzDecoder::new(File::open(&filename).unwrap()))
        };
    }
}

