
use jseqio::reader::*;
use jseqio::writer::*;
use std::str;

use std::env;

fn main(){
    let args: Vec<String> = env::args().collect();
    let infile = args[1].clone();
    let outdir = args[2].clone();

    let mut reader = DynamicFastXReader::from_file(&infile).unwrap();

    let mut prev: Vec<u8> = vec![];

    let mut writer: Option<DynamicFastXWriter> = None;

    while let Some(rec) = reader.read_next().unwrap(){
        let mut tokens = rec.head.split(|x| *x == b' ');
        let first = tokens.next().unwrap();
        let accession = first.split(|x| *x == b'.').next().unwrap();
        if accession != prev.as_slice(){
            // New accession number. Open a new fasta writer for that.

            let outfile = format!("{}/{}.fna", outdir, str::from_utf8(accession).unwrap());

            if let Some(mut w) = writer {
                w.flush().unwrap();
            }

            writer = Some(DynamicFastXWriter::new_to_file(&outfile).unwrap());
            eprintln!("Processing accession {}", String::from_utf8(accession.to_vec()).unwrap());
        }
        writer.as_mut().unwrap().write(&rec).unwrap();
        prev = accession.to_owned();
    }

}