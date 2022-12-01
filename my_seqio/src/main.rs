use std::io;
use std::fs::File;
use std::io::BufReader;
use std::fmt;
use std::str;

enum InputMode{
    FASTA,
    FASTQ,
}

struct FastXReader<R: io::BufRead>{
    inputmode: InputMode,
    input: R,
    seq_buf: Vec<u8>,
    head_buf: Vec<u8>,
    qual_buf: Vec<u8>,
    plus_buf: Vec<u8>, // For the fastq plus-line
}

#[derive(Debug)]
struct SeqRecord<'a>{
    seq: &'a [u8],
    head: &'a [u8],
    qual: Option<&'a [u8]>, // If FASTA, this is None
}

impl<'a> fmt::Display for SeqRecord<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
               "SeqRecord{{\n  Head: {}\n  Seq: {}\n  Qual: {}\n}}", 
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
    fn next(&mut self) -> Option<SeqRecord>{
        // Just single-line FASTQ for now
        self.seq_buf.clear();
        self.head_buf.clear();
        self.qual_buf.clear();
        self.plus_buf.clear();

        let bytes_read = self.input.read_until(b'\n', &mut self.head_buf); // Read header line
        if bytes_read.unwrap() == 0 {return None} // End of stream

        self.input.read_until(b'\n', &mut self.seq_buf); // Read sequence line
        self.input.read_until(b'\n', &mut self.plus_buf); // Read plus-line
        self.input.read_until(b'\n', &mut self.qual_buf); // Read the quality line


        return Some(SeqRecord{seq: self.seq_buf.split_last().unwrap().1, // split_last trims trailing newline
                              head: self.head_buf.split_last().unwrap().1, 
                              qual: Some(self.qual_buf.split_last().unwrap().1)});
    }

    fn new(input: R, mode: InputMode) -> Self{
        FastXReader{inputmode: mode,
                    input: input,
                    seq_buf: Vec::<u8>::new(),
                    head_buf: Vec::<u8>::new(),
                    qual_buf: Vec::<u8>::new(),
                    plus_buf: Vec::<u8>::new()}
    }

}

// Todo: Use the lending iterator crate


fn main() {
    let input = BufReader::new(File::open(&"reads.fastq").unwrap());
    let mut reader = FastXReader::new(input, InputMode::FASTA);
    loop{
        if let Some(record) = reader.next(){
            println!("{}", record);
        } else { break };
    }
}
