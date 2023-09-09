use seq_tools::*;
use jseqio::reader::*;
use jseqio::writer::DynamicFastXWriter;

use std::env;

fn main(){
    let args: Vec<String> = env::args().collect();
    let infile = args[1].clone();
    let out_prefix = args[2].clone();
    let n_subsamples: u64 = args[3].clone().parse().unwrap();

    let reader = DynamicFastXReader::from_file(&infile).unwrap();

    let n_seqs = count_sequences(reader); // Consumes the reader
    let mut prev_file = infile;

    for i in 0u64..n_subsamples {
        let howmany = n_seqs / 2u64.pow(i.try_into().unwrap()); // n / 2^i
        if howmany == 0 { break }
        let outfile = format!("{}-{}.fna.gz", out_prefix, howmany.to_string());

        eprintln!("Creating file {}", outfile);

        let mut reader = DynamicFastXReader::from_file(&prev_file).unwrap();
        let mut writer = DynamicFastXWriter::new_to_file(&outfile).unwrap();

        // copy `howmany` records to a new file.
        for _ in 0u64..howmany{
            let rec = reader.read_next().unwrap();
            writer.write(&rec.expect("Not enough records in file"));
        }

        prev_file = outfile.clone();

    }
}