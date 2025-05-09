use jseqio::{reader::*, record::*, writer::*, reverse_complement_in_place};

mod histogram;
pub mod trim_adapters;

use rand_chacha::rand_core::SeedableRng;
use std::cmp::{max,min};
use sha2::{Sha256, Digest};

struct LengthIterator<'a>{
    reader: &'a mut DynamicFastXReader,
}

impl<'a> Iterator for LengthIterator<'a>{
    type Item = i64;

    fn next(&mut self) -> Option<Self::Item>{
        let rec = self.reader.read_next().unwrap();
        rec.map(|r| r.seq.len() as i64)
    }
}

pub fn extract_region(mut reader: DynamicFastXReader, mut writer: DynamicFastXWriter, start: usize, end: usize){
    let rec = reader.read_next().unwrap().expect("First sequence not found");
    let seq = rec.seq;
    if end >= seq.len() {
        panic!("Ending point {} past the end of a sequence of length {}", end, seq.len());
    }

    let region_seq = &seq[start..end+1]; // End is inclusive
    let region_qual = rec.qual.map(|q| &q[start..end+1]);

    let region_rec = RefRecord{seq: region_seq, qual: region_qual, head: rec.head};
    writer.write_ref_record(&region_rec).unwrap();
}

pub fn extract_reads_by_names(reader: DynamicFastXReader, names: &Vec<String>){
    let filetype = reader.filetype();

    let db = reader.into_db().unwrap();

    // There may be multiple records with the same name, like in interleaved paired-end fastq data

    let mut name_to_seqs = std::collections::BTreeMap::<Vec<u8>, Vec<RefRecord>>::new();
    for name in names{
        name_to_seqs.insert(name.as_bytes().to_owned(), vec![]);
    }
    
    for rec in db.iter(){
        if name_to_seqs.contains_key(rec.name()){
            name_to_seqs.get_mut(rec.name()).unwrap().push(rec);
        }
    }

    // Print in order to stdout
    let mut writer = jseqio::writer::DynamicFastXWriter::new_to_stdout(filetype, jseqio::CompressionType::None);
    for name in names {
        for rec in name_to_seqs.get(name.as_bytes()).unwrap(){
            writer.write(rec).unwrap();
        }
    }

}

pub fn extract_reads_by_ranks(mut reader: DynamicFastXReader, ranks: &Vec<usize>){
    let filetype = reader.filetype();

    // Print in order to stdout
    let mut writer = jseqio::writer::DynamicFastXWriter::new_to_stdout(filetype, jseqio::CompressionType::None);

    if ranks.is_sorted() {
        // Can stream the records
        let mut seq_idx = 0_usize;
        let mut ranks_idx = 0_usize;
        while let Some(rec) = reader.read_next().unwrap() {
            while seq_idx == ranks[ranks_idx] { // While-loop so we are okay with duplicate ranks
                writer.write(&rec).unwrap();
                ranks_idx += 1;
                if ranks_idx == ranks.len() { // Done
                    return;
                }
            }
            seq_idx += 1;
        }
        panic!("Error: did not find read with rank {}", ranks[ranks_idx]);
    } else {
        // Can not stream the reads because we want to print the raads in the order
        // they come in the ranks vector.

        let db = reader.into_db().unwrap(); // Read all reads to memory
        for rank in ranks {
            let rec = db.get(*rank);
            writer.write(&rec).unwrap();
        }

    }

}

pub fn print_lengths(reader: &mut DynamicFastXReader){
    while let Some(rec) = reader.read_next().unwrap() {
        println!("{}", rec.seq.len());
    }
}

pub fn gc_content(reader: &mut DynamicFastXReader){
    let mut n_GC = 0_usize;
    let mut n_AT = 0_usize;
    let mut n_other = 0_usize;
    while let Some(rec) = reader.read_next().unwrap() {
        for c in rec.seq {
            match c.to_ascii_uppercase() {
                b'G' => n_GC += 1, 
                b'C' => n_GC += 1, 
                b'A' => n_AT += 1,
                b'T' => n_AT += 1,
                _ => n_other += 1
            }
        }
    }
    println!("{} nucleotides ignored ({}%)", n_other, n_other as f64 / (n_other + n_GC + n_AT) as f64 * 100.0); 
    println!("GC content: {:.2}%", n_GC as f64 / (n_GC + n_AT) as f64 * 100.0);

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

    while let Some(rec) = reader.read_next().unwrap() {
        total_length += rec.seq.len() as u64;
        number_of_sequences += 1;
        max_seq_len = max(max_seq_len, rec.seq.len() as u64);
        min_seq_len = min(min_seq_len, rec.seq.len() as u64);

        // Check quality values if they exist
        if let Some(qual) = rec.qual{
            for q in qual{
                let x = q - 0x21; // Fastq quality bytes start from 0x21
                min_quality_value = min(min_quality_value, x as u64);
                max_quality_value = max(max_quality_value, x as u64);
                sum_of_quality_values += x as u64;
            }
        }
    }

    println!("Number of nucleotides: {}", total_length);
    println!("Number of sequences: {}", number_of_sequences);
    println!("Maximum sequence length: {}", max_seq_len);
    println!("Minimum sequence length: {}", min_seq_len);
    println!("Average sequence length: {}", total_length as f64 / number_of_sequences as f64);
    if min_quality_value != u64::MAX{
        println!("Maximum quality value: {}", max_quality_value);
        println!("Minimum quality value: {}", min_quality_value);
        println!("Average quality value: {}", sum_of_quality_values as f64 / total_length as f64);
    }

}

// Removes sequenes that have exactly the same nucleotides. The headers need not match.
pub fn remove_duplicates(reader: &mut DynamicFastXReader, writer: &mut DynamicFastXWriter){
    let mut seen: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new(); // Hash values of seen sequences
    let mut hasher = Sha256::new();
    while let Some(rec) = reader.read_next().unwrap(){
        hasher.update(rec.seq);
        let hashvalue = hasher.finalize_reset();
        if !seen.contains(hashvalue.as_slice()){
            writer.write(&rec).unwrap();
            seen.insert(hashvalue.to_vec());
        }
    }
}


pub fn print_length_histogram(reader: &mut DynamicFastXReader, min: i64, max: i64, n_bins: i64){
    let it = LengthIterator{reader};
    histogram::print_histogram(it, min, max, n_bins);
}

pub fn count_sequences(mut input: DynamicFastXReader) -> u64{
    let mut count = 0u64;
    while input.read_next().unwrap().is_some(){
        count += 1;
    }
    count
}

// Returns a random permutation of [0..n_elements)
fn get_random_permutation(n_elements: usize, seed_option: Option<u64>) -> Vec<usize> {

    // Generate random seed from current time if not given
    let seed = match seed_option{
        Some(s) => s,
        None => {
            // Generate a random seed from the current time
            let now = std::time::SystemTime::now();
            let since_epoch = now.duration_since(std::time::UNIX_EPOCH).unwrap();
            (since_epoch.as_nanos() % 0xFFFFFFFFFFFFFFFF) as u64
        }
    };

    // Put the seed into a 32-byte array
    let mut seed_array: [u8; 32] = [0; 32];
    for i in 0..8{
        seed_array[i] = ((seed >> (i * 8)) & 0xFF) as u8;
    }

    // Create a random number generator with the fixed seed
    let mut rng = rand_chacha::ChaCha20Rng::from_seed(seed_array);

    eprintln!("Generating random numbers...");
    let mut v: Vec<(u64, usize)> = vec![]; // Pairs (random u64, index)

    for i in 0..n_elements{
        let r = rand::RngCore::next_u64(&mut rng);
        v.push((r, i));
    }

    eprintln!("Sorting random numbers...");
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    eprintln!("Sorting done");

    // Collect the permutation
    v.iter().map(|x| x.1).collect()
}

// Needs two input readers to the same data because needs
// to pass over the data twice. Seed is the random seed. If not given, a seed is generated from the current time.
pub fn random_subsample(input1: DynamicFastXReader, input2: DynamicFastXReader, out: &mut DynamicFastXWriter, fraction: f64, seed: Option<u64>, paired_interleaved: bool){
    let n_seqs = count_sequences(input1) as usize; // Consumes the input

    let mut subsample_seqs: usize = (n_seqs as f64 * fraction) as usize;
    if paired_interleaved && subsample_seqs % 2 > 0{
        subsample_seqs -= 1; // Get to an even number
    }

    random_subsample_howmany(input2, out, n_seqs, subsample_seqs, seed, paired_interleaved);
}

pub fn get_subsample_keep_marks(n_seqs: usize, subsample_seqs: usize, seed: Option<u64>) -> Vec<u8>{
    let perm = get_random_permutation(n_seqs, seed);

    let mut keep_marks: Vec<u8> = vec![0u8; perm.len()];
    for id in perm.iter().take(subsample_seqs){
        keep_marks[*id] = 1;
    }
    keep_marks
}

pub fn get_subsample_pair_keep_marks(n_seqs: usize, subsample_seqs: usize, seed: Option<u64>) -> Vec<u8>{
    if subsample_seqs % 2 != 0{
        panic!("Error: the number of sequences to subsample must be even when subsampling paired-end interleaved data");
    }
    if n_seqs % 2 != 0{
        panic!("Error: paired-end interleaved data has an odd number of sequences");
    }

    let n_total_pairs = n_seqs / 2;
    let n_sample_pairs = subsample_seqs / 2;

    let perm = get_random_permutation(n_total_pairs, seed);

    let mut keep_marks: Vec<u8> = vec![0u8; n_seqs];
    for id in perm.iter().take(n_sample_pairs){
        keep_marks[2 * id] = 1;
        keep_marks[2 * id + 1] = 1;
    }
    keep_marks
}


// Seed is the random seed. If not given, a seed is generated from the current time.
pub fn random_subsample_howmany(mut input: DynamicFastXReader, out: &mut DynamicFastXWriter, total_seqs: usize, subsample_seqs: usize, seed: Option<u64>, paired_interleaved: bool){
    if paired_interleaved && total_seqs % 2 != 0{
        panic!("Error: the number of sequences must be even when subsampling paired-end interleaved data");
    }

    let keep_marks = match paired_interleaved{
        false => get_subsample_keep_marks(total_seqs, subsample_seqs, seed),
        true => get_subsample_pair_keep_marks(total_seqs, subsample_seqs, seed),
    };

    let mut seq_idx = 0;
    while let Some(rec) = input.read_next().unwrap(){
        if keep_marks[seq_idx] == 1{
            out.write(&rec).unwrap();
        }
        seq_idx += 1;
    }

    eprintln!("Done");
}

pub fn convert(input: &mut DynamicFastXReader, output: &mut DynamicFastXWriter){
    let mut dummy_qual_values: Vec<u8> = vec![]; // A buffer for dummy quality values for fasta -> fastq conversion 
    while let Some(mut rec) = input.read_next().unwrap(){
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
        output.write(&rec).unwrap();
    }   
}

pub fn reverse_complement(input: &mut DynamicFastXReader, output: &mut DynamicFastXWriter){
    let mut dummy_qual_values: Vec<u8> = vec![]; // A buffer for dummy quality values for fasta -> fastq conversion 
    while let Some(mut rec) = input.read_next_mut().unwrap(){
        if matches!(rec.qual, None){
            // Potentially doing Fasta to Fastq conversion.
            // Put dummy quality values to rec.qual.
            while dummy_qual_values.len() < rec.seq.len(){
                dummy_qual_values.push(b'I');
                // 'I' is the maximum quality value from most sequencers.
                // Some software may break if they see quality values larger than 'I'.
                // Hence, we use 'I' as the a dummy value.
            }
            rec.qual = Some(&mut dummy_qual_values.as_mut_slice()[0..rec.seq.len()]);
        }
        reverse_complement_in_place(rec.seq);
        rec.qual.as_mut().map(|q| q.reverse()); // Also reverse the quality values

        output.write(&rec.into_shared_ref()).unwrap();
    }   
}

pub fn concatenate(input: &mut DynamicFastXReader, output: &mut DynamicFastXWriter, header: &[u8]){
    let mut seq_concat = Vec::<u8>::new();
    let mut qual_concat = Vec::<u8>::new();

    while let Some(rec) = input.read_next().unwrap(){
        if let Some(qual) = rec.qual{
            qual_concat.extend_from_slice(qual)
        }
        seq_concat.extend_from_slice(rec.seq);
    }

    let rec_out = jseqio::record::RefRecord{
        head: header, 
        seq: &seq_concat, 
        qual: if !qual_concat.is_empty() {Some(&qual_concat)} else {None}};

    output.write(&rec_out).unwrap();
}

pub fn trim(input: &mut DynamicFastXReader, output: &mut DynamicFastXWriter, from_start: usize, from_end: usize, min_final_len: usize){
    let mut n_deleted: u64 = 0;
    while let Some(mut rec) = input.read_next().unwrap(){
        if rec.seq.len() >= from_start + from_end + min_final_len{
            // Trimming leaves at least one nucleotide
            rec.seq = &rec.seq[from_start .. rec.seq.len() - from_end];
            if let Some(qual) = rec.qual{
                // Quality values are present -> trim those too
                rec.qual = Some(&qual[from_start .. qual.len() - from_end]);
            }
            output.write(&rec).unwrap();
        } else{
            // Delete this sequence
            n_deleted += 1;
        }
    }
    if n_deleted > 0 {
        eprintln!("Deleted {} sequences whose final length would have been below the minimum length {}", n_deleted, min_final_len);
    }
}


pub fn get_reader(args: &clap::ArgMatches) -> Result<DynamicFastXReader, Box<dyn std::error::Error>>{
    let filename = args.get_one::<String>("input");

    if let Some(infile) = filename {
        // From file
        DynamicFastXReader::from_file(&infile)
    } else {
        // From stdin
        DynamicFastXReader::from_stdin()
    }
}

pub fn get_writer(args: &clap::ArgMatches) -> DynamicFastXWriter{
    let filename = args.get_one::<String>("output");

    if let Some(outfile) = filename {
        // From file
        DynamicFastXWriter::new_to_file(&outfile).unwrap()
    } else {
        // To stdout
        let is_fasta = args.get_flag("fasta-out");
        let is_fastq = args.get_flag("fastq-out");
        let compression_type = match args.get_flag("gzip-out"){
            true => jseqio::CompressionType::Gzip,
            false => jseqio::CompressionType::None,
        };

        if !is_fasta && !is_fastq {
            panic!(
                "Error: must give --fasta-out or --fastq-out and possibly --gzip-out if writing to stdout."
            );
        };

        let filetype = if is_fastq {jseqio::FileType::FASTQ} else {jseqio::FileType::FASTA};
        DynamicFastXWriter::new_to_stdout(filetype, compression_type)
    }
}

mod tests{
    use super::*;

    #[test]
    fn unit_test_paired_subsample(){

        let n_seqs = 100;
        let subsample_seqs = 50;
        let seed = Some(1234);

        let marks = get_subsample_pair_keep_marks(n_seqs, subsample_seqs, seed);
        let marks_sum = marks.iter().fold(0, |sum, x| sum + *x as usize);

        assert_eq!(marks_sum, subsample_seqs);

        // Check that sampling is pairwise
        for i in 0..n_seqs{
            if i % 2 == 0 && marks[i] == 1{
                assert_eq!(marks[i+1], 1);
            }
        }
    }
}