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

pub struct TemplateTest<T: std::io::Read>{
    reader: seq_io::fastq::Reader<T>
}

impl SeqStream for GzippedFastqStream{

    fn read_all(&mut self){
        let mut sum: i64 = 0;
        while let Some(record) = self.file_gz_reader.next() {
            let record = record.expect("Error reading record");
            sum += record.seq().len() as i64;
        }
        println!("{}", sum);
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
        let mut sum: i64 = 0;
        while let Some(record) = self.file_reader.next() {
            let record = record.expect("Error reading record");
            sum += record.seq().len() as i64;
        }
        println!("{}", sum);
    }
}

impl RegularFastqStream{
    pub fn new(filename: &String) -> RegularFastqStream{
        return RegularFastqStream{
            file_reader: Reader::new(File::open(&filename).unwrap())
        };
    }
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

