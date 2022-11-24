extern crate flate2;

use crate::fastq_streams::{FastaStream, FastqStream, SeqStream};

mod fastq_streams;

use clap::{Arg, ArgAction, Command};
use fastq_streams::MyRecord;
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{self, Write};

struct SeqReader {
    stream: Box<dyn SeqStream>,
}

impl SeqReader {
    // New from file
    pub fn new(filename: &String) -> SeqReader {
        if filename.ends_with("fastq.gz") {
            return SeqReader {
                stream: Box::new(FastqStream::<GzDecoder<File>>::new_from_file(filename)),
            };
        } else if filename.ends_with("fastq") {
            return SeqReader {
                stream: Box::new(FastqStream::<File>::new_from_file(filename)),
            };
        } else if filename.ends_with("fna.gz") {
            return SeqReader {
                stream: Box::new(FastaStream::<GzDecoder<File>>::new_from_file(filename)),
            };
        } else if filename.ends_with("fna") {
            return SeqReader {
                stream: Box::new(FastaStream::<File>::new_from_file(filename)),
            };
        } else {
            panic!("Could not determine the format of file {}", filename);
        }
    }

    // New from stdin
    pub fn new_from_stdin(fastq: bool, gzipped: bool) -> SeqReader {
        if fastq && gzipped {
            return SeqReader {
                stream: Box::new(FastqStream::<GzDecoder<std::io::Stdin>>::new(
                    GzDecoder::new(io::stdin()),
                )),
            };
        } else if fastq && !gzipped {
            return SeqReader {
                stream: Box::new(FastqStream::<std::io::Stdin>::new(io::stdin())),
            };
        } else if !fastq && gzipped {
            return SeqReader {
                stream: Box::new(FastaStream::<GzDecoder<std::io::Stdin>>::new(
                    GzDecoder::new(io::stdin()),
                )),
            };
        } else if !fastq && !gzipped {
            return SeqReader {
                stream: Box::new(FastaStream::<std::io::Stdin>::new(io::stdin())),
            };
        } else {
            panic!("This line should never be reached");
        }
    }

    // Returns None if no more records
    pub fn read_next(&mut self) -> Option<MyRecord>{
        return self.stream.next_record()
    }

}

fn print_all_to_stdout(reader: &mut SeqReader){
    loop{
        let next = reader.read_next();
        match next {
            Some(rec) => {
                std::io::stdout().write_all(b"Header: ").ok();
                std::io::stdout().write_all(rec.header.as_slice()).ok();
                std::io::stdout().write_all(b"\n").ok();
                std::io::stdout().write_all(b"Sequence: ").ok();
                std::io::stdout().write_all(rec.seq.as_slice()).ok();
                std::io::stdout().write_all(b"\n").ok();
                std::io::stdout().write_all(b"Quality: ").ok();
                std::io::stdout().write_all(match rec.qual {Some(x) => x, None => b"Not available".to_vec()}.as_slice()).ok();
                std::io::stdout().write_all(b"\n").ok();
            }
            None => break
        }
    }
}

fn print_stats(reader: &mut SeqReader){
    let mut total_length: usize = 0; // i32 to test what happens with overflow
    loop{
        match reader.read_next() {
            Some(rec) => total_length += rec.seq.len(),
            None => break
        }
    }
    println!("Number of nucleotides: {}", total_length);
}

fn main() {
    let matches = Command::new("Fasta/fastq parsing")
        .version("0.1.0")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .about("Fasta/fastq parsing")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .help("Input filename"),
        )
        .arg(
            Arg::new("fasta")
                .short('a')
                .long("fasta")
                .action(ArgAction::SetTrue)
                .help("Parse in fasta format"),
        )
        .arg(
            Arg::new("fastq")
                .short('q')
                .long("fastq")
                .action(ArgAction::SetTrue)
                .help("Parse in fastq format"),
        )
        .arg(
            Arg::new("gzip")
                .short('g')
                .long("gzip")
                .action(ArgAction::SetTrue)
                .help("Read gzipped input"),
        ).arg(
            Arg::new("stats")
                .short('s')
                .long("stats")
                .action(ArgAction::SetTrue)
                .help("Print stats about the input."),
        )
        .get_matches();

    if let Some(infile) = matches.get_one::<String>("input") {
        let mut reader = SeqReader::new(&infile);
        print_all_to_stdout(&mut reader);
    } else {
        let is_fasta = matches.get_flag("fasta");
        let is_fastq = matches.get_flag("fastq");
        let is_gzip = matches.get_flag("gzip");
        if is_fasta && is_fastq {
            println!("Error: can't give both fasta and fastq flags.");
            std::process::exit(-1)
        }
        if !is_fasta && !is_fastq {
            println!(
                "Error: must give --fasta or --fastq and possibly --gzip if reading from stdin."
            );
            std::process::exit(-1)
        }
        let mut reader = SeqReader::new_from_stdin(is_fastq, is_gzip);
        if matches.get_flag("stats") {
            print_stats(&mut reader);
        }
    }
}
