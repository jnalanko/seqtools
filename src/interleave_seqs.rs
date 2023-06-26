use seq_tools::*;
use my_seqio::reader::DynamicFastXReader;
use my_seqio::writer::{DynamicFastXWriter, self};
use std::collections::HashSet;
use std::str;
use std::cmp::max;
use std::fs::File;
use std::io::{BufWriter};

use std::env;

fn main(){
    let args: Vec<String> = env::args().collect();
    let seqs1 = args[1].clone();
    let seqs2 = args[2].clone();
    let outfile = args[3].clone();

    let mut reader1 = DynamicFastXReader::new_from_file(&seqs1);
    let mut reader2 = DynamicFastXReader::new_from_file(&seqs2);

    let mut writer = DynamicFastXWriter::new_to_file(&outfile);

    while let Some(rec1) = reader1.read_next(){
        let rec2 = reader2.read_next().expect("File 2 has fewer records than file 1.");
        writer.write(rec1);
        writer.write(rec2);
    }

    writer.flush();
    if let Some(_) = reader2.read_next(){
        panic!("File 2 has more records than file 1.")
    }

}