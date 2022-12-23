use my_seqio::reader::DynamicFastXReader;
use my_seqio::writer::DynamicFastXWriter;

mod histogram;

use std::io::{Write,BufWriter};
use rand::Rng;
use std::cmp::{max,min};

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

pub fn print_stats(reader: &mut DynamicFastXReader){
    let mut total_length: u64 = 0;
    let mut number_of_sequences: u64 = 0;
    let mut max_seq_len: u64 = 0;
    let mut min_seq_len: u64 = u64::MAX;

    // Quality value statistics, if exist
    let mut max_quality_value: u64 = 0;
    let mut min_quality_value: u64 = u64::MAX;
    let mut sum_of_quality_values: u64 = 0;

    loop{
        match reader.read_next() {
            Some(rec) => {
                total_length += rec.seq.len() as u64;
                number_of_sequences += 1;
                max_seq_len = max(max_seq_len, rec.seq.len() as u64);
                min_seq_len = min(min_seq_len, rec.seq.len() as u64);

                // Check quality values is they exist
                if let Some(qual) = rec.qual{
                    for q in qual{
                        let x = q - 0x21; // Fastq quality bytes start from 0x21
                        min_quality_value = min(min_quality_value, x as u64);
                        max_quality_value = max(max_quality_value, x as u64);
                        sum_of_quality_values += x as u64;
                    }
                }
            },
            None => break
        }
    }

    println!("Number of nucleotides: {}", total_length);
    println!("Number of sequences: {}", number_of_sequences);
    println!("Maximum sequence length: {}", max_seq_len);
    println!("Minimum sequence length: {}", min_seq_len);
    println!("Average sequence length: {}", total_length as f64 / number_of_sequences as f64);
    if max_quality_value != u64::MAX{
        println!("Maximum quality value: {}", max_quality_value);
        println!("Minimum quality value: {}", min_quality_value);
        println!("Average quality value: {}", sum_of_quality_values as f64 / total_length as f64);
    }

}


pub fn print_length_histogram(reader: &mut DynamicFastXReader, min: i64, max: i64, n_bins: i64){
    let it = LengthIterator{reader: reader};
    histogram::print_histogram(it, min, max, n_bins);
}

pub fn count_sequences(mut input: DynamicFastXReader) -> u64{
    let mut count = 0u64;
    while let Some(_) = input.read_next(){
        count += 1;
    }
    return count;
}

// Returns a random permutation of [0..n_elements)
fn get_random_permutation(n_elements: usize) -> Vec<usize> {
    let mut v: Vec<(f64, usize)> = vec![]; // Random number from 0 to 1, index
    let mut rng = rand::thread_rng();

    for i in 0..n_elements{
        let r = rng.gen_range(0.0..1.0);
        v.push((r, i));
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Build the permutation
    let mut p: Vec<usize> = vec![];
    for i in 0..n_elements{
        p.push(v[i].1);
    }
    p
}

// Needs two input reader to the same data because needs
// to pass over the data twice.
pub fn random_subsample(input1: DynamicFastXReader, input2: DynamicFastXReader, out: &mut DynamicFastXWriter, fraction: f64){
    let n_seqs = count_sequences(input1) as usize; // Consumes the input
    let subsample_seqs: usize = (n_seqs as f64 * fraction) as usize;

    random_subsample_howmany(input2, out, n_seqs, subsample_seqs);
}

pub fn random_subsample_howmany(mut input: DynamicFastXReader, out: &mut DynamicFastXWriter, total_seqs: usize, sumsample_seqs: usize){
    let perm = get_random_permutation(total_seqs); 

    let mut keep_marks: Vec<u8> = vec![0u8; perm.len()];
    for id in perm.iter().take(sumsample_seqs){
        keep_marks[*id as usize] = 1;
    }

    let mut seq_idx = 0;
    while let Some(rec) = input.read_next(){
        if keep_marks[seq_idx] == 1{
            out.write(rec);
        }
        seq_idx += 1;
    }

    eprintln!("Done");
}

pub fn convert(input: &mut DynamicFastXReader, output: &mut DynamicFastXWriter){
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

pub fn trim(input: &mut DynamicFastXReader, output: &mut DynamicFastXWriter, from_start: usize, from_end: usize){
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
    if n_deleted > 0 {
        eprintln!("Deleted {} sequences as too short to trim", n_deleted);
    }
}


pub fn get_reader(args: &clap::ArgMatches) -> DynamicFastXReader{
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
                "Error: must give --fasta-in or --fastq-in and possibly --gzip-in if reading from stdin."
            );
        };
        let filetype = if is_fastq {my_seqio::FileType::FASTQ} else {my_seqio::FileType::FASTA};
        DynamicFastXReader::new_from_stdin(filetype, is_gzip)
    }
}

pub fn get_writer(args: &clap::ArgMatches) -> DynamicFastXWriter{
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