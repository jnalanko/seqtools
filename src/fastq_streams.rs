use flate2::read::GzDecoder;
use seq_io::fasta::Reader as seqio_fasta_reader;
use seq_io::fasta::Record as seqio_fasta_record;
use seq_io::fastq::Reader as seqio_fastq_reader;
use seq_io::fastq::Record as seqio_fastq_record;
use std::fs::File;

pub struct MyRecord{
    pub seq: Vec<u8>,
    pub header: Vec<u8>,
    pub qual: Option<Vec<u8>>
}

// Fasta or fastq stream
pub trait SeqStream {
    fn read_all(&mut self);
    fn next_record(&mut self) -> Option<MyRecord>;
}

//
//
// FASTQ START
//
//

// Template class for a fastq stream that takes the raw data stream as a template parameter
pub struct FastqStream<T: std::io::Read> {
    reader: seqio_fastq_reader<T>,
}

impl<T: std::io::Read> SeqStream for FastqStream<T> {
    fn read_all(&mut self) {
        let mut sum: i64 = 0;
        while let Some(record) = self.reader.next() {
            let record = record.expect("Error reading record");
            sum += record.seq().len() as i64;
        }
        println!("{}", sum);
    }

    fn next_record(&mut self) -> Option<MyRecord>{
        let opt = self.reader.next();
        return match opt{
            Some(res) => match res {
                Ok(rec) => Some(MyRecord{seq: rec.seq().to_vec(), 
                                         header: rec.head().to_vec(), 
                                         qual: Some(rec.qual().to_vec())}),
                Err(e) => panic!("Error reading fastq file: {}", e.to_string())
            },
            None => None
        };
    }
}

impl FastqStream<File> {
    pub fn new_from_file(filename: &String) -> FastqStream<File> {
        return FastqStream::<File> {
            reader: seqio_fastq_reader::new(File::open(&filename).unwrap()),
        };
    }
}

impl FastqStream<GzDecoder<File>> {
    pub fn new_from_file(filename: &String) -> FastqStream<GzDecoder<File>> {
        return FastqStream::<GzDecoder<File>> {
            reader: seqio_fastq_reader::new(GzDecoder::new(File::open(&filename).unwrap())),
        };
    }
}

impl<T: std::io::Read> FastqStream<T> {
    pub fn new(stream: T) -> FastqStream<T> {
        return FastqStream::<T> {
            reader: seqio_fastq_reader::<T>::new(stream),
        };
    }
}

//
//
// FASTQ START
//
//

impl<T: std::io::Read> FastaStream<T> {
    pub fn new(stream: T) -> FastaStream<T> {
        return FastaStream::<T> {
            reader: seqio_fasta_reader::<T>::new(stream),
        };
    }
}

impl<T: std::io::Read> SeqStream for FastaStream<T> {
    fn read_all(&mut self) {
        let mut sum: i64 = 0;
        while let Some(record) = self.reader.next() {
            let record = record.expect("Error reading record");
            sum += record.seq().len() as i64;
        }
        println!("{}", sum);
    }

    fn next_record(&mut self) -> Option<MyRecord>{
        let opt = self.reader.next();
        return match opt{
            Some(res) => match res {
                Ok(rec) => Some(MyRecord{seq: rec.seq().to_vec(), 
                                         header: rec.head().to_vec(), 
                                         qual: None}),
                Err(e) => panic!("Error reading fastq file: {}", e.to_string())
            },
            None => None
        };
    }
}

impl FastaStream<File> {
    pub fn new_from_file(filename: &String) -> FastaStream<File> {
        return FastaStream::<File> {
            reader: seqio_fasta_reader::new(File::open(&filename).unwrap()),
        };
    }
}

impl FastaStream<GzDecoder<File>> {
    pub fn new_from_file(filename: &String) -> FastaStream<GzDecoder<File>> {
        return FastaStream::<GzDecoder<File>> {
            reader: seqio_fasta_reader::new(GzDecoder::new(File::open(&filename).unwrap())),
        };
    }
}


// Template class for a fasta stream that takes the raw data stream as a template parameter
pub struct FastaStream<T: std::io::Read> {
    reader: seqio_fasta_reader<T>,
}