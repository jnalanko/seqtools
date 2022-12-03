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
    head: &'a [u8],    
    seq: &'a [u8],
    qual: Option<&'a [u8]>, // If FASTA, this is None
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
    fn next(&mut self) -> Option<SeqRecord>{
        // Just single-line FASTQ for now
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

#[cfg(test)]
mod tests {
    use super::*;

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
            "#################################################",
            "#################################################",
            "#################################################"
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
}


fn main() {
    let input = BufReader::new(File::open(&"reads_trunc.fastq").unwrap());
    let mut reader = FastXReader::new(input, InputMode::FASTA);
    loop{
        if let Some(record) = reader.next(){
            println!("{}", record);
        } else { break };
    }
}
