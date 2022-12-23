use seq_tools::*;
use my_seqio::reader::DynamicFastXReader;
use my_seqio::writer::DynamicFastXWriter;

use std::env;

fn main(){
    let args: Vec<String> = env::args().collect();
    let infile = args[1].clone();
    let out_prefix = args[2].clone();
    let n_subsamples: u64 = args[3].clone().parse().unwrap();

    let reader = DynamicFastXReader::new_from_file(&infile);

    let n_seqs = count_sequences(reader); // Consumes the reader
    let mut prev_file = infile.clone();

    for i in 0u64..n_subsamples {
        let howmany = n_seqs / 2u64.pow(i.try_into().unwrap()); // n / 2^i
        if howmany == 0 { break }
        let outfile = format!("{}-{}.fna.gz", out_prefix, howmany.to_string());

        let mut reader = DynamicFastXReader::new_from_file(&prev_file);
        let mut writer = DynamicFastXWriter::new_to_file(&outfile);

        // copy `howmany` records to a new file.
        for _ in 0u64..howmany{
            let rec = reader.read_next();
            writer.write(rec.unwrap());
        }

        prev_file = outfile.clone();

    }
}