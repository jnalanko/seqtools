extern crate flate2;

mod histogram;

use my_seqio::reader::DynamicFastXReader;
use my_seqio::writer::DynamicFastXWriter;

use std::io::{Write,BufWriter};
use rand::Rng;
use std::cmp::{max,min};

mod cli;

struct LengthIterator<'a>{
    reader: &'a mut DynamicFastXReader,
}

impl<'a> Iterator for LengthIterator<'a>{
    type Item = i64;

    fn next(&mut self) -> Option<Self::Item>{
        let rec = self.reader.read_next();
        match rec{
            None => None,
            Some(r) => Some(r.seq.len() as i64),
        }
    }
}

fn print_stats(reader: &mut DynamicFastXReader){
    let mut total_length: u64 = 0;
    let mut number_of_sequences: u64 = 0;
    let mut max_seq_len: u64 = 0;
    let mut min_seq_len: u64 = 1e20 as u64;
    loop{
        match reader.read_next() {
            Some(rec) => {
                total_length += rec.seq.len() as u64;
                number_of_sequences += 1;
                max_seq_len = max(max_seq_len, rec.seq.len() as u64);
                min_seq_len = min(min_seq_len, rec.seq.len() as u64);
            },
            None => break
        }
    }
    println!("Number of nucleotides: {}", total_length);
    println!("Number of sequences: {}", number_of_sequences);
    println!("Maximum sequence length: {}", max_seq_len);
    println!("Minimum sequence length: {}", min_seq_len);
    println!("Average sequence length: {}", total_length as f64 / number_of_sequences as f64);
}

fn print_length_histogram(reader: &mut DynamicFastXReader, min: i64, max: i64, n_bins: i64){
    let it = LengthIterator{reader: reader};
    histogram::print_histogram(it, min, max, n_bins);
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

    eprintln!("Done");
}

fn convert(input: &mut DynamicFastXReader, output: &mut DynamicFastXWriter){
    let mut dummy_qual_values: Vec<u8> = vec![]; // A buffer for dummy quality values for fasta -> fastq conversion 
    while let Some(mut rec) = input.read_next(){
        if matches!(rec.qual, None){
            // Potentially doing Fasta to Fastq conversion.
            // Put dummy quality values to rec.qual.
            while dummy_qual_values.len() < rec.seq.len(){
                dummy_qual_values.push(b'I');
                // 'I' is the maximum quality value from most sequencers.
                // Some software may break if they see quality values larger than 'I'.
                // Hence, we use 'I' as the a dummy value.
            }
            rec.qual = Some(&dummy_qual_values.as_slice()[0..rec.seq.len()]);
        }
        output.write(rec);
    }   
}

fn trim(input: &mut DynamicFastXReader, output: &mut DynamicFastXWriter, from_start: usize, from_end: usize){
    let mut n_deleted: u64 = 0;
    while let Some(mut rec) = input.read_next(){
        if rec.seq.len() > from_start + from_end{
            // Trimming leaves at least one nucleotide
            rec.seq = &rec.seq[from_start .. rec.seq.len() - from_end];
            if let Some(qual) = rec.qual{
                // Quality values are present -> trim those too
                rec.qual = Some(&qual[from_start .. qual.len() - from_end]);
            }
            output.write(rec);
        } else{
            // Delete this sequence
            n_deleted += 1;
        }
    }
    eprintln!("Deleted {} sequences as too short to trim", n_deleted);
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
        Some(("convert", _)) => { 
            let mut reader = get_reader(&matches);
            let mut writer = get_writer(&matches);
            convert(&mut reader, &mut writer);
        }
        Some(("trim", sub_matches)) => { 
            let mut reader = get_reader(&matches);
            let mut writer = get_writer(&matches);
            let from_start: usize = sub_matches.get_one::<String>("from-start").unwrap().parse().unwrap();
            let from_end: usize = sub_matches.get_one::<String>("from-end").unwrap().parse().unwrap();
            trim(&mut reader, &mut writer, from_start, from_end);
        }
        _ => {}
    };
}
