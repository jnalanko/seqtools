
use jseqio::reader::DynamicFastXReader;
use jseqio::writer::{DynamicFastXWriter};






use std::env;

fn main(){
    let args: Vec<String> = env::args().collect();
    let seqs1 = args[1].clone();
    let seqs2 = args[2].clone();
    let outfile = args[3].clone();

    let mut reader1 = DynamicFastXReader::new_from_file(&seqs1).unwrap();
    let mut reader2 = DynamicFastXReader::new_from_file(&seqs2).unwrap();

    let mut writer = DynamicFastXWriter::new_to_file(&outfile).unwrap();

    while let Some(rec1) = reader1.read_next().unwrap(){
        let rec2 = reader2.read_next().expect("File 2 has fewer records than file 1.").unwrap();
        writer.write(&rec1);
        writer.write(&rec2);
    }

    writer.flush();
    if let Some(_) = reader2.read_next().unwrap(){
        panic!("File 2 has more records than file 1.")
    }

}