
use std::io;
use std::io::BufWriter;
use std::io::Write;
use std::fs::File;
use flate2::Compression;
use flate2::write::GzEncoder;

use crate::FileType;
use crate::record::SeqRecord;
use crate::record::Record;
use crate::figure_out_file_format;

pub struct FastXWriter<W: Write>{
    pub filetype: FileType,
    pub output: BufWriter<W>,
}

pub trait SeqRecordWriter{
    // Can not take a Record trait object because then we can't
    // for some reason put a SeqRecordWriter into a box.
    // So we take the header, sequence and quality values as slices.
    fn write(&mut self, head: &[u8], seq: &[u8], qual: Option<&[u8]>);
}

pub struct DynamicFastXWriter {
    stream: Box<dyn SeqRecordWriter>,
}

impl DynamicFastXWriter{

    pub fn write<Rec: Record>(&mut self, rec: Rec){
        self.stream.write(rec.head(), rec.seq(), rec.qual());
    }

    pub fn new_to_stream<W: Write + 'static>(stream: W, filetype: FileType) -> Self{
        let writer = FastXWriter::<W>::new(stream, filetype);
        DynamicFastXWriter {stream: Box::new(writer)}
    }

    // Write to a file
    pub fn new_to_file(filename: &String) -> Self {
        let output = File::create(&filename).unwrap();
        match figure_out_file_format(&filename.as_str()){
            (FileType::FASTQ, true) =>{
                let gzencoder = GzEncoder::<File>::new(output, Compression::fast());
                Self::new_to_stream(gzencoder, FileType::FASTQ)
            },
            (FileType::FASTQ, false) => {
                Self::new_to_stream(output, FileType::FASTQ)
            },
            (FileType::FASTA, true) => {
                let gzencoder = GzEncoder::<File>::new(output, Compression::fast());
                Self::new_to_stream(gzencoder, FileType::FASTA)
            },
            (FileType::FASTA, false) => {
                Self::new_to_stream(output, FileType::FASTA)
            },
        }
    }

    pub fn new_to_stdout(filetype: FileType, gzipped: bool) -> Self {
        if gzipped {
            Self::new_to_stream(GzEncoder::new(io::stdout(), Compression::fast()), filetype)
        } else {
            Self::new_to_stream(io::stdout(), filetype)
        }
    }
}

impl<W: Write> FastXWriter<W>{
    pub fn write<Rec: Record>(&mut self, rec: &Rec){
        match &self.filetype{
            FileType::FASTA => {
                self.output.write(b">").expect("Error writing output");
                self.output.write(rec.head()).expect("Error writing output");
                self.output.write(b"\n").expect("Error writing output");
                self.output.write(rec.seq()).expect("Error writing output");
                self.output.write(b"\n").expect("Error writing output");
            }
            FileType::FASTQ => {
                self.output.write(b"@").expect("Error writing output");
                self.output.write(rec.head()).expect("Error writing output");
                self.output.write(b"\n").expect("Error writing output");
                self.output.write(rec.seq()).expect("Error writing output");
                self.output.write(b"\n+\n").expect("Error writing output");
                self.output.write(rec.qual().expect("Quality values missing")).expect("Error writing output");
                self.output.write(b"\n").expect("Error writing output");
            }
        }
    }

    pub fn new(output: W, filetype: FileType) -> Self{
        Self{
            filetype: filetype,
            output: BufWriter::<W>::new(output)
        }
    }

    pub fn flush(&mut self){
        self.output.flush().expect("Error flushing output stream");
    }
}

impl<W: Write> SeqRecordWriter for FastXWriter<W>{
    fn write(&mut self, head: &[u8], seq: &[u8], qual: Option<&[u8]>){
        let rec = SeqRecord{head, seq, qual};
        self.write(&rec);
    }
}