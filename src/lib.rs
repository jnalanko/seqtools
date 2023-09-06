use jseqio::{reader::DynamicFastXReader, record::RefRecord, record::OwnedRecord, record::Record};
use jseqio::writer::DynamicFastXWriter;

mod histogram;
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

pub fn extract_reads_by_names(reader: DynamicFastXReader, names: &Vec<String>){
    let filetype = reader.filetype();
    let mut names_hashset = std::collections::HashSet::new();
    for name in names{
        names_hashset.insert(name.as_bytes());
    }
    
    let db = reader.into_db().unwrap();
    let mut found: Vec<RefRecord> = db.iter().filter(|rec| names_hashset.contains(rec.name())).collect();

    // Sort by the order in the given ranks
    let mut names_order = std::collections::HashMap::<&[u8], usize>::new();
    for (i, name) in names.iter().enumerate(){
        names_order.insert(name.as_bytes(), i);
    }

    found.sort_by_key(|x| names_order.get(x.name()).unwrap());

    // Print in order to stdout
    let mut writer = jseqio::writer::DynamicFastXWriter::new_to_stdout(filetype, false);
    for rec in found.iter(){
        writer.write(rec);
    }

}

pub fn extract_reads_by_ranks(mut reader: DynamicFastXReader, ranks: &Vec<usize>){
    let mut ranks_hashset = std::collections::HashSet::new();
    for r in ranks{
        ranks_hashset.insert(r);
    }

    let mut found = Vec::<(usize, OwnedRecord)>::new(); // Key, record
    let mut current_rank = 0_usize;
    while let Some(rec) = reader.read_next().unwrap() {
        if ranks_hashset.contains(&current_rank) {
            found.push((current_rank, rec.to_owned()));
        }
        current_rank += 1;
    }

    if found.len() != ranks_hashset.len(){
        panic!("Error: Could not find all reads");
    }

    // Sort by the order in the given ranks
    let mut ranks_order = std::collections::HashMap::<usize, usize>::new();
    for (i, r) in ranks.iter().enumerate(){
        ranks_order.insert(*r, i);
    }

    found.sort_by_key(|x| ranks_order.get(&x.0).unwrap());

    // Print in order to stdout
    let mut writer = jseqio::writer::DynamicFastXWriter::new_to_stdout(reader.filetype(), false);
    for (_, rec) in found.iter(){
        writer.write(rec);
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
            writer.write(&rec);
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
pub fn random_subsample(input1: DynamicFastXReader, input2: DynamicFastXReader, out: &mut DynamicFastXWriter, fraction: f64, seed: Option<u64>){
    let n_seqs = count_sequences(input1) as usize; // Consumes the input
    let subsample_seqs: usize = (n_seqs as f64 * fraction) as usize;

    random_subsample_howmany(input2, out, n_seqs, subsample_seqs, seed);
}


// Seed is the random seed. If not given, a seed is generated from the current time.
pub fn random_subsample_howmany(mut input: DynamicFastXReader, out: &mut DynamicFastXWriter, total_seqs: usize, subsample_seqs: usize, seed: Option<u64>){
    let perm = get_random_permutation(total_seqs, seed); 

    let mut keep_marks: Vec<u8> = vec![0u8; perm.len()];
    for id in perm.iter().take(subsample_seqs){
        keep_marks[*id] = 1;
    }

    let mut seq_idx = 0;
    while let Some(rec) = input.read_next().unwrap(){
        if keep_marks[seq_idx] == 1{
            out.write(&rec);
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
        output.write(&rec);
    }   
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
            output.write(&rec);
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
        let is_gzip = args.get_flag("gzip-out");
        if is_fasta && is_fastq {
            panic!("Error: can't give both fasta and fastq flags.");
        }
        if !is_fasta && !is_fastq {
            panic!(
                "Error: must give --fasta-out or --fastq-out and possibly --gzip-out if writing to stdout."
            );
        };
        let filetype = if is_fastq {jseqio::FileType::FASTQ} else {jseqio::FileType::FASTA};
        DynamicFastXWriter::new_to_stdout(filetype, is_gzip)
    }
}