extern crate flate2;

use crate::fastq_streams::{FastaStream, FastqStream, SeqStream};

mod fastq_streams;

use clap::{Arg, ArgAction, Command};
use flate2::read::GzDecoder;
use std::env;
use std::ffi::OsStr;
use std::fs::File;
use std::io;
use std::path::Path;


fn main() {
    let args: Vec<String> = env::args().collect();
    //let mut m_args = args.clone();
    let seq = args[1].clone(); // Take a shallow copy of the string
    let mut bytes = seq.into_bytes(); // Turn the string into a byte vector (consumes string)
    let n = bytes.len();

    // Complement
    for i in 0..n {
        bytes[i] = match bytes[i]{
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => {
                println!("Warning: non-ACGT character: '{}'", bytes[i] as char); 
                bytes[i]
            }
        };
    }

    // Reverse
    for i in 0..(n/2){
        let temp = bytes[i];
        bytes[i] = bytes[n-1-i];
        bytes[n-1-i] = temp;
    }

    println!("{}", std::str::from_utf8(bytes.as_slice()).expect(""));
}
