use std::io;
use std::fs::File;

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

struct SeqRecord<'a>{
    seq: &'a [u8],
    head: &'a [u8],
    qual: Option<&'a [u8]>, // If FASTA, this is None
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

        return Some(SeqRecord{seq: self.seq_buf.as_slice(),
                              head: self.head_buf.as_slice(), 
                              qual: Some(self.qual_buf.as_slice())});

    }
}

// Todo: Use the lending iterator crate


fn main() {
    let f = File::open(&"reads.fastq").unwrap();
}
