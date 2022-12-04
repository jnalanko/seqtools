
use std::io;
use std::io::BufReader;
use std::io::BufRead;
use std::io::BufWriter;
use std::io::Write;
use std::fmt;
use std::str;
use std::fs::File;
use flate2::read::MultiGzDecoder;

pub enum InputMode{
    FASTA,
    FASTQ,
}
pub enum OutputMode{
    FASTA,
    FASTQ,
}

pub struct FastXReader<R: io::Read>{
    inputmode: InputMode,
    input: BufReader<R>,
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

impl<R: io::Read> FastXReader<R>{
    pub fn next(&mut self) -> Option<SeqRecord>{
        if matches!(self.inputmode, InputMode::FASTQ){
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
                panic!("FASTQ quality line has different length than sequence line")
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

    pub fn new(input: R, mode: InputMode) -> Self{
        FastXReader{inputmode: mode,
                    input: BufReader::new(input),
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
}

// Implement common SeqStream trait for all
// FastXReaders over the generic parameter R.
impl<R: io::Read> SeqRecordProducer for FastXReader<R>{
    fn next(&mut self) -> Option<SeqRecord>{
        self.next()
    }
}

pub struct DynamicFastXReader {
    stream: Box<dyn SeqRecordProducer>,
}

// A class that contains a dynamic trait object for different
// types of input streams.
impl DynamicFastXReader {

    // Need to constrain + 'static because boxed things always need to have a static
    // lifetime.
    pub fn new_from_input_stream<R: io::Read + 'static>(r: R, mode: InputMode) -> Self{
        let reader = FastXReader::<R>::new(r, mode);
        DynamicFastXReader {stream: Box::new(reader)}
    }

    // New from file
    pub fn new_from_file(filename: &String) -> Self {
        let input = File::open(&filename).unwrap();
        if filename.ends_with("fastq.gz") {
            let gzdecoder = MultiGzDecoder::<File>::new(input);
            Self::new_from_input_stream(gzdecoder, InputMode::FASTQ)
        } else if filename.ends_with("fastq") {
            Self::new_from_input_stream(input, InputMode::FASTQ)
        } else if filename.ends_with("fna.gz") {
            let gzdecoder = MultiGzDecoder::<File>::new(input);
            Self::new_from_input_stream(gzdecoder, InputMode::FASTA)
        } else if filename.ends_with("fna") {
            Self::new_from_input_stream(input, InputMode::FASTA)
        } else {
            panic!("Could not determine the format of file {}", filename);
        }
    }

    // New from stdin
    pub fn new_from_stdin(fastq: bool, gzipped: bool) -> Self {
        let mode = if fastq {InputMode::FASTQ} else {InputMode::FASTA};
        if gzipped {
            Self::new_from_input_stream(MultiGzDecoder::new(io::stdin()), mode)
        } else {
            Self::new_from_input_stream(io::stdin(), mode)
        }
    }

    // Returns None if no more records
    pub fn read_next(&mut self) -> Option<SeqRecord>{
        return self.stream.next()
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
        let mut reader = FastXReader::new(input, InputMode::FASTQ);

        let mut seqs_read = 0;
        loop{
            if let Some(record) = reader.next(){
                assert_eq!(record.head, headers[seqs_read].as_bytes());
                assert_eq!(record.seq, seqs[seqs_read].as_bytes());
                assert_eq!(record.qual.unwrap(), quals[seqs_read].as_bytes());
                seqs_read += 1;
            } else { break };
        }
        assert_eq!(seqs_read, n_seqs);
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
        let mut reader = FastXReader::new(input, InputMode::FASTA);

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
        let mut writer = FastXWriter::<Vec<u8>>::new(out_buf, OutputMode::FASTA);

        for rec in owned_records.iter() {
            writer.write(rec);
        }

        writer.flush();
        let written_data = writer.output.into_inner().unwrap();

        // This written data may not exactly equal the original data,
        // because the length of FASTA sequence lines is not fixed.
        // Read the records back from written data and compare to originals.

        let mut reader2 = FastXReader::new(written_data.as_slice(), InputMode::FASTA);
        let mut seqs_read2 = 0;
        loop{
            if let Some(record) = reader2.next(){
                dbg!(&record);
                assert_eq!(record.head, headers[seqs_read].as_bytes());
                assert_eq!(record.seq, seqs[seqs_read].as_bytes());
                assert_eq!(record.qual, None);
                seqs_read2 += 1;
            } else { break };
        }
        assert_eq!(seqs_read2, n_seqs);

    }
}

pub struct FastXWriter<W: Write>{
    outputmode: OutputMode,
    output: BufWriter<W>,
}

impl<W: Write> FastXWriter<W>{
    pub fn write<Rec: Record>(&mut self, rec: &Rec){
        match &self.outputmode{
            OutputMode::FASTA => {
                self.output.write(b">").expect("Error writing output");
                self.output.write(rec.head()).expect("Error writing output");
                self.output.write(b"\n").expect("Error writing output");
                self.output.write(rec.seq()).expect("Error writing output");
                self.output.write(b"\n").expect("Error writing output");
            }
            OutputMode::FASTQ => {
                self.output.write(b"@").expect("Error writing output");
                self.output.write(rec.head()).expect("Error writing output");
                self.output.write(b"\n").expect("Error writing output");
                self.output.write(rec.seq()).expect("Error writing output");
                self.output.write(b"\n+\n").expect("Error writing output");
                self.output.write(rec.qual().expect("Quality values missing")).expect("Error writing output");
            }
        }
    }

    pub fn new(output: W, mode: OutputMode) -> Self{
        Self{
            outputmode: mode,
            output: BufWriter::<W>::new(output)
        }
    }

    pub fn flush(&mut self){
        self.output.flush().expect("Error flushing output stream");
    }
}