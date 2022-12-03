extern crate flate2;

use my_seqio;

use clap::{Arg, ArgAction, Command};
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{self, Write};
use std::env;

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

/*
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
*/

fn print_stats(reader: &mut SeqReader){
    let mut total_length: usize = 0;
    let mut number_of_sequences: usize = 0;
    loop{
        match reader.read_next() {
            Some(rec) => {
                total_length += rec.seq.len();
                number_of_sequences += 1;  
            },
            None => break
        }
    }
    println!("Number of nucleotides: {}", total_length);
    println!("Number of sequences: {}", number_of_sequences);
}

fn print_length_histogram(reader: &mut SeqReader, min: i64, max: i64, n_bins: i64){

    let mut counters: Vec<i64> = vec![0; n_bins as usize];
    let bin_width = (max-min+1) / n_bins;
    loop{
        match reader.read_next() {
            Some(rec) => {
                let len = rec.seq.len() as i64;
                let mut bin = (len - min as i64) / bin_width;

                // Clamp to [0, n_bins-1]
                bin = std::cmp::max(0, bin);
                bin = std::cmp::min(n_bins-1, bin);

                counters[bin as usize] += 1;
            },
            None => break
        }
    }

    let max_counter: i64 = *counters.iter().max().unwrap();
    let n_columns: i64 = 40;

    for (i, c) in counters.iter().enumerate(){
        let n_chars = ((*c as f64 / max_counter as f64) * n_columns as f64) as i64;
        print!("{}\t", i*(bin_width as usize));
        std::io::stdout().write_all(vec![b'#'; n_chars as usize].as_slice()).ok();
        println!();
    }
}

enum ReaderInput{
    FromFile{filename: String},
    FromStdIn{is_fastq: bool, is_gzipped: bool} // Is fasta if not fastq
}

fn get_reader(input: ReaderInput) -> SeqReader{
    match input{
        ReaderInput::FromFile{filename} => SeqReader::new(&filename),
        ReaderInput::FromStdIn{is_fastq, is_gzipped} => SeqReader::new_from_stdin(is_fastq, is_gzipped)
    }
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
                .help("Input filename")
                .global(true),
        )
        .arg(
            Arg::new("fasta")
                .short('a')
                .long("fasta")
                .action(ArgAction::SetTrue)
                .help("Parse in fasta format")
                .global(true),
        )
        .arg(
            Arg::new("fastq")
                .short('q')
                .long("fastq")
                .action(ArgAction::SetTrue)
                .help("Parse in fastq format")
                .global(true),
        )
        .arg(
            Arg::new("gzip")
                .short('g')
                .long("gzip")
                .action(ArgAction::SetTrue)
                .help("Read gzipped input")
                .global(true),
        ).subcommand(Command::new("length-histogram")
            .about("Print the length histogram of the sequences.")
            .arg(
                Arg::new("min")
                    .default_value("0")
                    .help("Minimum value")
            ).arg(
                Arg::new("max")
                    .default_value("1000")
                    .help("Maximum value")
            ).arg(
                Arg::new("nbins")
                    .default_value("20")
                    .help("Number of bins")
            )
        ).subcommand(Command::new("stats")
            .about("Print stats about the input.")
        )
        .get_matches();


    // Flag consistency check
    let mut reader = if let Some(infile) = matches.get_one::<String>("input") {
        get_reader(ReaderInput::FromFile{filename: infile.clone()})
    } else {
        // From stdin
        let is_fasta = matches.get_flag("fasta");
        let is_fastq = matches.get_flag("fastq");
        let is_gzip = matches.get_flag("gzip");
        if is_fasta && is_fastq {
            println!("Error: can't give both fasta and fastq flags.");
            std::process::exit(-1);
        }
        if !is_fasta && !is_fastq {
            println!(
                "Error: must give --fasta or --fastq and possibly --gzip if reading from stdin."
            );
            std::process::exit(-1);
        };
        get_reader(ReaderInput::FromStdIn{is_fastq: is_fastq, is_gzipped: is_gzip})
    };

    match matches.subcommand() {
        Some(("length-histogram", sub_matches)) => { 
            let min: i64 = sub_matches.get_one::<String>("min").unwrap().parse::<i64>().unwrap();
            let max: i64 = sub_matches.get_one::<String>("max").unwrap().parse::<i64>().unwrap();
            let nbins: i64 = sub_matches.get_one::<String>("nbins").unwrap().parse::<i64>().unwrap();
            print_length_histogram(&mut reader, min as i64, max as i64, nbins as i64);
        }
        Some(("stats", _)) => { 
            print_stats(&mut reader);
        }
        _ => {}
    };
}
