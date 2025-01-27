
use jseqio::reader::*;
use jseqio::writer::*;
use std::env;

fn main(){
    let args: Vec<String> = env::args().collect();
    let infile = args[1].clone();
    let out1 = args[2].clone();
    let out2 = args[3].clone();

    let mut reader = DynamicFastXReader::from_file(&infile).unwrap();

    let mut writer1 = DynamicFastXWriter::new_to_file(&out1).unwrap();
    let mut writer2 = DynamicFastXWriter::new_to_file(&out2).unwrap();

    let mut n_seqs_read = 0_usize;
    let mut flag = true;
    while let Some(rec) = reader.read_next().unwrap(){
        if flag {
            writer1.write_ref_record(&rec).unwrap();
        } else {
            writer2.write_ref_record(&rec).unwrap();
        }
        flag = !flag;
        n_seqs_read += 1;
    }

    if n_seqs_read % 2 != 0 {
        panic!("Odd number of sequences in input file.")
    }

}