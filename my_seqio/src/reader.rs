use std::io;
use std::fs::File;
use std::io::BufReader;
use flate2::read::MultiGzDecoder;
use crate::FileType;
use crate::record::SeqRecord;
use crate::figure_out_file_format;

pub struct FastXReader<R: io::BufRead>{
    pub filetype: FileType,
    pub input: R,
    pub seq_buf: Vec<u8>,
    pub head_buf: Vec<u8>,
    pub qual_buf: Vec<u8>,
    pub plus_buf: Vec<u8>, // For the fastq plus-line
    pub fasta_temp_buf: Vec<u8>, // Stores the fasta header read in the previous iteration
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

pub struct DynamicFastXReader {
    stream: Box<dyn SeqRecordProducer>,
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