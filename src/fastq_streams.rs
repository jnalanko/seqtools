use std::io;
use flate2::read::GzDecoder;
use seq_io::fastq::Reader as seqio_fastq_reader;
use seq_io::fastq::Record as seqio_fastq_record;
use seq_io::fasta::Reader as seqio_fasta_reader;
use seq_io::fasta::Record as seqio_fasta_record;
use std::fs::File;

// Fasta or fastq stream
pub trait SeqStream{
    fn read_all(&mut self);
}

// Template class for a fastq stream that takes the raw data stream as a template parameter
pub struct FastqStream<T: std::io::Read>{
    reader: seqio_fastq_reader<T>
}

// Template class for a fasta stream that takes the raw data stream as a template parameter
pub struct FastaStream<T: std::io::Read>{
    reader: seqio_fasta_reader<T>
}

impl<T: std::io::Read> SeqStream for FastqStream<T>{

    fn read_all(&mut self){
        let mut sum: i64 = 0;
        while let Some(record) = self.reader.next() {
            let record = record.expect("Error reading record");
            sum += record.seq().len() as i64;
        }
        println!("{}", sum);
    }
}

impl FastqStream<File>{
    pub fn new(filename: &String) -> FastqStream<File>{
        return FastqStream::<File>{
            reader: seqio_fastq_reader::new(File::open(&filename).unwrap())
        };
    }
}

impl FastqStream<std::io::Stdin>{
    pub fn new() -> FastqStream<std::io::Stdin>{
        return FastqStream::<std::io::Stdin>{
            reader: seqio_fastq_reader::new(io::stdin())
        };
    }
}

impl FastqStream<GzDecoder<File>>{
    pub fn new(filename: &String) -> FastqStream<GzDecoder<File>>{
        return FastqStream::<GzDecoder<File>>{
            reader: seqio_fastq_reader::new(GzDecoder::new(File::open(&filename).unwrap()))
        };
    }
}

impl<T: std::io::Read> SeqStream for FastaStream<T>{

    fn read_all(&mut self){
        let mut sum: i64 = 0;
        while let Some(record) = self.reader.next() {
            let record = record.expect("Error reading record");
            sum += record.seq().len() as i64;
        }
        println!("{}", sum);
    }
}

impl FastaStream<File>{
    pub fn new(filename: &String) -> FastaStream<File>{
        return FastaStream::<File>{
            reader: seqio_fasta_reader::new(File::open(&filename).unwrap())
        };
    }
}

impl FastaStream<GzDecoder<File>>{
    pub fn new(filename: &String) -> FastaStream<GzDecoder<File>>{
        return FastaStream::<GzDecoder<File>>{
            reader: seqio_fasta_reader::new(GzDecoder::new(File::open(&filename).unwrap()))
        };
    }
}


