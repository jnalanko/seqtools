extern crate flate2;

use my_seqio::{DynamicFastXReader};
use std::io::{Write};

mod cli;

fn print_stats(reader: &mut DynamicFastXReader){
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

fn print_length_histogram(reader: &mut DynamicFastXReader, min: i64, max: i64, n_bins: i64){

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
        print!("{}\t", (min + (i as i64)*bin_width) as usize);
        std::io::stdout().write_all(vec![b'#'; n_chars as usize].as_slice()).ok();
        println!();
    }
}

enum ReaderInput{
    FromFile{filename: String},
    FromStdIn{is_fastq: bool, is_gzipped: bool} // Is fasta if not fastq
}

fn get_reader(input: ReaderInput) -> DynamicFastXReader{
    match input{
        ReaderInput::FromFile{filename} => DynamicFastXReader::new_from_file(&filename),
        ReaderInput::FromStdIn{is_fastq, is_gzipped} => DynamicFastXReader::new_from_stdin(is_fastq, is_gzipped)
    }
}

fn main() {

    let matches = cli::build_cli().get_matches();

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
