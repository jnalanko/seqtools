use seq_tools::*;
use jseqio::reader::DynamicFastXReader;
use jseqio::writer::DynamicFastXWriter;
use std::collections::HashSet;
use std::str;
use std::cmp::max;

use std::env;

fn hash_kmers(filename: &String, k: usize) -> HashSet<Vec<u8>>{
    let mut reader = DynamicFastXReader::new_from_file(&filename).unwrap();
    let mut kmers: HashSet<Vec<u8>> = HashSet::new();
    while let Some(rec) = reader.read_next().unwrap(){
        let seq = rec.seq;
        let n = seq.len() as i64;
        let k = k as i64;
        let m = max(0, n-k+1) as usize; // Number of kmers in the sequence
        for i in 0..m{
            let kmer = &seq[i..i + k as usize];
            kmers.insert(kmer.to_vec());
        }
    }
    kmers
}

fn main(){
    let args: Vec<String> = env::args().collect();
    let seqs1 = args[1].clone();
    let seqs2 = args[2].clone();
    let k: usize = args[3].parse().unwrap();

    let K1 = hash_kmers(&seqs1, k);
    let K2 = hash_kmers(&seqs2, k);

    let intersection_size = K1.intersection(&K2).count();
    let union_size = K1.union(&K2).count();

    eprintln!("File 1 has {} distinct {}-mers", K1.len(), k);
    eprintln!("File 2 has {} distinct {}-mers", K2.len(), k);
    eprintln!("The intersection has {} distinct {}-mers", intersection_size, k);
    eprintln!("The union has {} distinct {}-mers", union_size, k);
    eprintln!("Jaccard index: {}", intersection_size as f64 / union_size as f64);

}