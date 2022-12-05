
use std::io;
use std::io::BufReader;
use std::io::BufRead;
use std::io::BufWriter;
use std::io::Write;
use std::fmt;
use std::str;
use std::fs::File;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder; // todo: flate2::bufwrite?

#[derive(Copy, Clone)]
pub enum FileType{
    FASTA,
    FASTQ,
}

pub struct FastXReader<R: io::BufRead>{
    filetype: FileType,
    input: R,
    seq_buf: Vec<u8>,
    head_buf: Vec<u8>,
    qual_buf: Vec<u8>,
    plus_buf: Vec<u8>, // For the fastq plus-line
    fasta_temp_buf: Vec<u8>, // Stores the fasta header read in the previous iteration
}

pub trait Record{
    fn head(&self) -> &[u8];
    fn seq(&self) -> &[u8];
    fn qual(&self) -> Option<&[u8]>;
}

#[derive(Debug)]
pub struct SeqRecord<'a>{
    pub head: &'a [u8],    
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>, // If FASTA, this is None
}

#[derive(Debug)]
pub struct OwnedSeqRecord{
    pub head: Vec<u8>,    
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>, // If FASTA, this is None
}

impl<'a> Record for SeqRecord<'a>{
    fn head(&self) -> &[u8]{self.head}
    fn seq(&self) -> &[u8]{self.seq}
    fn qual(&self) -> Option<&[u8]>{self.qual}
}

impl<'a> Record for OwnedSeqRecord{
    fn head(&self) -> &[u8]{self.head.as_slice()}
    fn seq(&self) -> &[u8]{self.seq.as_slice()}
    fn qual(&self) -> Option<&[u8]>{
        match &self.qual{
            Some(q) => return Some(q.as_slice()),
            None => None,
        }
    }
}

impl<'a> SeqRecord<'a>{
    pub fn to_owned(&self) -> OwnedSeqRecord{
        OwnedSeqRecord { 
            head: self.head.to_vec(), 
            seq: self.seq.to_vec(), 
            qual: match self.qual {
                Some(q) => Some(q.to_vec()), 
                None => None
            }
        }
    }
}

impl<'a> fmt::Display for SeqRecord<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
               "SeqRecord{{ \n  Head: {}\n  Seq:  {}\n  Qual: {}\n}}", 
               str::from_utf8(self.head).unwrap(),
               str::from_utf8(self.seq).unwrap(),
               match self.qual{
                   Some(q) => str::from_utf8(q).unwrap(),
                   None => "", // No quality values
               }
               
        )
    }
}

impl<R: io::BufRead> FastXReader<R>{
    pub fn next(&mut self) -> Option<SeqRecord>{
        if matches!(self.filetype, FileType::FASTQ){
            // FASTQ format

            self.seq_buf.clear();
            self.head_buf.clear();
            self.qual_buf.clear();
            self.plus_buf.clear();

            // Read header line
            let bytes_read = self.input.read_until(b'\n', &mut self.head_buf);
            if bytes_read.expect("I/O error.") == 0 {return None} // End of stream

            // Read sequence line
            let bytes_read = self.input.read_until(b'\n', &mut self.seq_buf);
            if bytes_read.expect("I/O error.") == 0 {
                panic!("FASTQ sequence line missing."); // File can't end here
            }
            
            // read +-line
            let bytes_read = self.input.read_until(b'\n', &mut self.plus_buf);
            if bytes_read.expect("I/O error.") == 0 {
                panic!("FASTQ + line missing."); // File can't end here
            }

            // read qual-line
            let bytes_read = self.input.read_until(b'\n', &mut self.qual_buf);
            let bytes_read = bytes_read.expect("I/O error.");
            if bytes_read == 0{ // File can't end here
                panic!("FASTQ quality line missing."); 
            } else if bytes_read != self.seq_buf.len(){
                panic!("FASTQ quality line has different length than sequence line ({} vs {})", bytes_read, self.seq_buf.len())
            }

            return Some(SeqRecord{head: self.head_buf.as_slice().strip_prefix(b"@").unwrap().strip_suffix(b"\n").unwrap(), 
                                seq: self.seq_buf.as_slice().strip_suffix(b"\n").unwrap(),
                                qual: Some(self.qual_buf.as_slice().strip_suffix(b"\n").unwrap())})
        }
        else{
            // FASTA format
            self.seq_buf.clear();
            self.head_buf.clear();

            // Read header line
            if self.fasta_temp_buf.len() == 0 {
                // This is the first record -> read header from input
                let bytes_read = self.input.read_until(b'\n', &mut self.head_buf);
                if bytes_read.expect("I/O error.") == 0 {return None} // End of stream
            } else{
                // Take stashed header from previous iteration
                self.head_buf.append(&mut self.fasta_temp_buf); // Also clears the temp buf
            }

            // Read sequence line
            loop{
                let bytes_read = self.input.read_until(b'\n', &mut self.fasta_temp_buf);
                match bytes_read{
                    Err(e) => panic!("{}",e), // File can't end here
                    Ok(bytes_read) => {
                        if bytes_read == 0{
                            // No more bytes left to read
                            if self.seq_buf.len() == 0{
                                // Stream ends with an empty sequence
                                panic!("Empty sequence in FASTA file");
                            }
                            break; // Ok, last record of the file
                        }

                        // Check if we read the header of the next read
                        let start = self.fasta_temp_buf.len() as isize - bytes_read as isize;
                        if self.fasta_temp_buf[start as usize] == b'>'{
                            // Found a header. Leave it to the buffer for the next iteration.
                            break;
                        } else{
                            // Found more sequence -> Append to self.seq_buf
                            self.seq_buf.append(&mut self.fasta_temp_buf); // Also clears the temp buf
                            self.seq_buf.pop(); // Trim newline (TODO: what if there is none?)
                        }
                    }
                }
            }

            return Some(SeqRecord{head: self.head_buf.as_slice().strip_prefix(b">").unwrap().strip_suffix(b"\n").unwrap(), 
                                seq: self.seq_buf.as_slice(), // Newlines are already trimmed before
                                qual: None});
        }
    }

    pub fn new(input: R, filetype: FileType) -> Self{
        FastXReader{filetype: filetype,
                    input: input,
                    seq_buf: Vec::<u8>::new(),
                    head_buf: Vec::<u8>::new(),
                    qual_buf: Vec::<u8>::new(),
                    plus_buf: Vec::<u8>::new(),
                    fasta_temp_buf: Vec::<u8>::new(),}
    }

}

// Trait for a stream returning SeqRecord objects.
pub trait SeqRecordProducer {
    fn next(&mut self) -> Option<SeqRecord>;
    fn filetype(&self )-> FileType; 
}

// Implement common SeqStream trait for all
// FastXReaders over the generic parameter R.
impl<R: io::BufRead> SeqRecordProducer for FastXReader<R>{
    fn next(&mut self) -> Option<SeqRecord>{
        self.next()
    }

    fn filetype(&self)-> FileType{
        self.filetype
    } 

}

pub struct DynamicFastXReader {
    stream: Box<dyn SeqRecordProducer>,
}

// Returns (file type, is_gzipped)
fn figure_out_file_format(filename: &str) -> (FileType, bool){
    let is_gzipped = filename.ends_with(".gz");
    let filename = if is_gzipped{
        &filename[0 .. filename.len()-3] // Drop the .gz suffix
    } else {&filename};
    let fasta_extensions = vec![".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa"];
    let fastq_extensions = vec![".fastq", ".fq"];
    if fasta_extensions.iter().any(|&suffix| filename.ends_with(suffix)){
        return (FileType::FASTA, is_gzipped);
    } else if fastq_extensions.iter().any(|&suffix| filename.ends_with(suffix)){
        return (FileType::FASTQ, is_gzipped);
    } else{
        panic!("Unkown file extension: {}", filename);
    }
}

// A class that contains a dynamic trait object for different
// types of input streams.
impl DynamicFastXReader {

    // Need to constrain + 'static because boxed things always need to have a static
    // lifetime.
    pub fn new_from_input_stream<R: io::BufRead + 'static>(r: R, filetype: FileType) -> Self{
        let reader = FastXReader::<R>::new(r, filetype);
        DynamicFastXReader {stream: Box::new(reader)}
    }

    // New from file
    pub fn new_from_file(filename: &String) -> Self {
        let input = File::open(&filename).unwrap();
        let (fileformat, gzipped) = figure_out_file_format(&filename.as_str());
        if gzipped{
            let gzdecoder = MultiGzDecoder::<File>::new(input);
            Self::new_from_input_stream(BufReader::new(gzdecoder), fileformat)
        } else{
            Self::new_from_input_stream(BufReader::new(input), fileformat)
        }
    }

    // New from stdin
    pub fn new_from_stdin(filetype: FileType, gzipped: bool) -> Self {
        if gzipped {
            Self::new_from_input_stream(BufReader::new(MultiGzDecoder::new(io::stdin())), filetype)
        } else {
            Self::new_from_input_stream(BufReader::new(io::stdin()), filetype)
        }
    }

    // Returns None if no more records
    pub fn read_next(&mut self) -> Option<SeqRecord>{
        return self.stream.next()
    }

    pub fn filetype(&self)-> FileType{
        self.stream.filetype()
    } 

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{cmp::min, process::Output};

    #[test]
    fn fastq() {
        let headers = vec!(
            "SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1",
            "SRR403017.2 HWUSI-EAS108E_0007:3:1:10327:976/1",
            "SRR403017.3 HWUSI-EAS108E_0007:3:1:13569:972/1");
        let seqs = vec!(
            "TTGGACCGGCGCAAGACGGACCAGNGCGAAAGCATTTGCCAAGAANNNN",
            "CAACTTTCTATCTGGCATTCCCTGNGGAGGAAATAGAATGCGCGCNNNN",
            "GATCGGAAGAGCACACGTCTGAACNCCAGTCACTTAGGCATCTCGNNNN",
        );
        let quals = vec!(
            "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQ",
            "RSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~####",
            "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        );

        let n_seqs = headers.len();
        let mut fastq_data: String = "".to_owned();
        for i in 0..n_seqs{
            fastq_data.push_str(format!("@{}\n", headers[i]).as_str());
            fastq_data.push_str(format!("{}\n", seqs[i]).as_str());
            fastq_data.push_str("+\n");
            fastq_data.push_str(format!("{}\n", quals[i]).as_str());
        }

        let input = BufReader::new(fastq_data.as_bytes());
        let mut reader = FastXReader::new(input, FileType::FASTQ);

        let mut owned_records: Vec<OwnedSeqRecord> = vec![];
        let mut seqs_read = 0;
        loop{
            if let Some(record) = reader.next(){
                assert_eq!(record.head, headers[seqs_read].as_bytes());
                assert_eq!(record.seq, seqs[seqs_read].as_bytes());
                assert_eq!(record.qual.unwrap(), quals[seqs_read].as_bytes());
                owned_records.push(record.to_owned());
                seqs_read += 1;
            } else { break };
        }
        assert_eq!(seqs_read, n_seqs);

        // Test writer
        let out_buf: Vec<u8> = vec![];
        let mut writer = FastXWriter::<Vec<u8>>::new(out_buf, FileType::FASTQ);

        for rec in owned_records.iter() {
            writer.write(rec);
        }

        writer.flush();
        let written_data = writer.output.into_inner().unwrap();

        // Read the records back from written data and compare to originals.

        let mut reader2 = FastXReader::new(written_data.as_slice(), FileType::FASTQ);
        let mut seqs_read2 = 0;
        loop{
            if let Some(record) = reader2.next(){
                dbg!(&record);
                assert_eq!(record.head, headers[seqs_read2].as_bytes());
                assert_eq!(record.seq, seqs[seqs_read2].as_bytes());
                assert_eq!(record.qual.unwrap(), quals[seqs_read2].as_bytes());
                seqs_read2 += 1;
            } else { break };
        }
        assert_eq!(seqs_read2, n_seqs);
    }

    #[test]
    fn fasta() {
        let headers: Vec<String> = vec!(
            "SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1".to_owned(),
            "SRR403017.2 HWUSI-EAS108E_0007:3:1:10327:976/1".to_owned(),
            "SRR403017.3 HWUSI-EAS108E_0007:3:1:13569:972/1".to_owned());
        let seqs: Vec<String> = vec!(
            "TTGGACCGGCGCAAGACGGACCAGNGCGAAAGCATTTGCCAAGAANNNN".to_owned(),
            "CAACTTTCTATCTGGCATTCCCTGNGGAGGAAATAGAATGCGCGCNNNN".to_owned(),
            "GATCGGAAGAGCACACGTCTGAACNCCAGTCACTTAGGCATCTCGNNNN".to_owned(),
        );

        fn split_seq_to_lines(seq: &String, line_length: usize) -> Vec<String>{
            let mut i: usize = 0;
            let mut lines = Vec::<String>::new();
            while line_length*i < seq.len(){
                lines.push(seq[line_length*i .. min(line_length*(i+1), seq.len())].to_owned());
                i += 1;
            }
            lines
        }

        let n_seqs = headers.len();
        let mut fasta_data: String = "".to_owned();
        for i in 0..n_seqs{
            fasta_data.push_str(format!(">{}\n", headers[i].as_str()).as_str());
            for line in split_seq_to_lines(&seqs[i], 11){
                // Line length 11 should make it so that the last line has a different
                // length than the other lines.
                fasta_data.push_str(format!("{}\n", line.as_str()).as_str());
            }
        }

        dbg!(&fasta_data);

        let input = BufReader::new(fasta_data.as_bytes());
        let mut reader = FastXReader::new(input, FileType::FASTA);

        let mut owned_records: Vec<OwnedSeqRecord> = vec![];
        let mut seqs_read = 0;
        loop{
            if let Some(record) = reader.next(){
                dbg!(&record);
                assert_eq!(record.head, headers[seqs_read].as_bytes());
                assert_eq!(record.seq, seqs[seqs_read].as_bytes());
                assert_eq!(record.qual, None);
                owned_records.push(record.to_owned());
                seqs_read += 1;
            } else { break };
        }
        assert_eq!(seqs_read, n_seqs);

        // Test writer
        let out_buf: Vec<u8> = vec![];
        let mut writer = FastXWriter::<Vec<u8>>::new(out_buf, FileType::FASTA);

        for rec in owned_records.iter() {
            writer.write(rec);
        }

        writer.flush();
        let written_data = writer.output.into_inner().unwrap();

        // This written data may not exactly equal the original data,
        // because the length of FASTA sequence lines is not fixed.
        // Read the records back from written data and compare to originals.

        let mut reader2 = FastXReader::new(written_data.as_slice(), FileType::FASTA);
        let mut seqs_read2 = 0;
        loop{
            if let Some(record) = reader2.next(){
                dbg!(&record);
                assert_eq!(record.head, headers[seqs_read2].as_bytes());
                assert_eq!(record.seq, seqs[seqs_read2].as_bytes());
                assert_eq!(record.qual, None);
                seqs_read2 += 1;
            } else { break };
        }
        assert_eq!(seqs_read2, n_seqs);

    }

    #[test]
    fn test_figure_out_file_format(){
        assert!(match figure_out_file_format("aa.fna") {(FileType::FASTA,false) => true, _ => false});
        assert!(match figure_out_file_format("aa.fq") {(FileType::FASTQ,false) => true, _ => false});
        assert!(match figure_out_file_format("bbb.fna.gz") {(FileType::FASTA,true) => true, _ => false});
        assert!(match figure_out_file_format("cc.fna.gz") {(FileType::FASTA,true) => true, _ => false});
        assert!(match figure_out_file_format(".fna.gz") {(FileType::FASTA,true) => true, _ => false});
        assert!(match figure_out_file_format(".fasta") {(FileType::FASTA,false) => true, _ => false});
        assert!(match figure_out_file_format(".fq") {(FileType::FASTQ,false) => true, _ => false});
    }
}

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