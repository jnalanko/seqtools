extern crate flate2;

use my_seqio::{DynamicFastXReader, FastXReader, DynamicFastXWriter, FastXWriter};
use std::io::{Write,BufWriter};
use rand::Rng;

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

// Needs two input reader to the same data because needs
// to pass over the data twice.
fn random_subsample(input1: &mut DynamicFastXReader, input2: &mut DynamicFastXReader, out: &mut DynamicFastXWriter, fraction: f64){
    let mut v: Vec<(f64, usize)> = vec![]; // Random number from 0 to 1, seq id
    let mut rng = rand::thread_rng();
    let mut seq_idx = 0;

    eprintln!("Assigning random numbers to sequences...");
    while let Some(_) = input1.read_next(){
        let r = rng.gen_range(0.0..1.0);
        v.push((r, seq_idx));
        seq_idx += 1;
    }

    eprintln!("{} sequences found", seq_idx);

    let howmany: usize = (v.len() as f64 * fraction) as usize;
    eprintln!("Subsampling {}% ({} sequences...)", fraction*100.0, howmany);

    v.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut keep_marks: Vec<u8> = vec![0u8; v.len()];
    for (_, id) in v.iter().take(howmany){
        keep_marks[*id as usize] = 1;
    }

    let mut seq_idx = 0;
    while let Some(rec) = input2.read_next(){
        if keep_marks[seq_idx] == 1{
            out.write(rec);
        }
        seq_idx += 1;
    }

    eprintln!("Done...");
}

fn get_reader(args: &clap::ArgMatches) -> DynamicFastXReader{
    let filename = args.get_one::<String>("input");

    if let Some(infile) = filename {
        // From file
        DynamicFastXReader::new_from_file(&infile)
    } else {
        // From stdin
        let is_fasta = args.get_flag("fasta-in");
        let is_fastq = args.get_flag("fastq-in");
        let is_gzip = args.get_flag("gzip-in");
        if is_fasta && is_fastq {
            panic!("Error: can't give both fasta and fastq flags.");
        }
        if !is_fasta && !is_fastq {
            panic!(
                "Error: must give --fasta-in or --fastq-in and possibly --gzip if reading from stdin."
            );
        };
        let filetype = if is_fastq {my_seqio::FileType::FASTQ} else {my_seqio::FileType::FASTA};
        DynamicFastXReader::new_from_stdin(filetype, is_gzip)
    }
}

fn get_writer(args: &clap::ArgMatches) -> DynamicFastXWriter{
    let filename = args.get_one::<String>("output");

    if let Some(outfile) = filename {
        // From file
        DynamicFastXWriter::new_to_file(&outfile)
    } else {
        // To stdout
        let is_fasta = args.get_flag("fasta-out");
        let is_fastq = args.get_flag("fastq-out");
        let is_gzip = args.get_flag("gzip-out");
        if is_fasta && is_fastq {
            panic!("Error: can't give both fasta and fastq flags.");
        }
        if !is_fasta && !is_fastq {
            panic!(
                "Error: must give --fasta-out or --fastq-out and possibly --gzip-out if writing to stdout."
            );
        };
        let filetype = if is_fastq {my_seqio::FileType::FASTQ} else {my_seqio::FileType::FASTA};
        DynamicFastXWriter::new_to_stdout(filetype, is_gzip)
    }
}


fn main() {

    let matches = cli::build_cli().get_matches();

    match matches.subcommand() {
        Some(("length-histogram", sub_matches)) => { 
            let mut reader = get_reader(&matches);
            let min: i64 = sub_matches.get_one::<String>("min").unwrap().parse::<i64>().unwrap();
            let max: i64 = sub_matches.get_one::<String>("max").unwrap().parse::<i64>().unwrap();
            let nbins: i64 = sub_matches.get_one::<String>("nbins").unwrap().parse::<i64>().unwrap();
            print_length_histogram(&mut reader, min as i64, max as i64, nbins as i64);
        }
        Some(("stats", _)) => { 
            let mut reader = get_reader(&matches);
            print_stats(&mut reader);
        }
        Some(("subsample", sub_matches)) => {
            if matches.get_one::<String>("input") == None {
                panic!("Can not subsample from stdin because we need to pass over the data twice.");
            }

            // Get two readers for two passes over the data
            let mut input1 = get_reader(&matches);
            let mut input2 = get_reader(&matches);
            let mut output = get_writer(&matches);
            let frac: f64 = sub_matches.get_one::<String>("fraction")
                .unwrap().parse::<f64>().unwrap();
            random_subsample(&mut input1, &mut input2, &mut output, frac);
        }
        _ => {}
    };
}
