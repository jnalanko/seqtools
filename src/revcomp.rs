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
    let seq = args[1].clone();
    let mut bytes = seq.into_bytes();
    let n = bytes.len();

    dbg!(n);

    // Reverse
    for i in 0..(n/2){
        let temp = bytes[i];
        bytes[i] = bytes[n-1-i];
        bytes[n-1-i] = temp;
    }

    //dbg!(bytes);

    // Complement
    for i in 0..n {
        if bytes[i] == (b'A'){
            bytes[i] = b'T';
        }
        else if bytes[i] == (b'C'){
            bytes[i] = b'G';
        }
        else if bytes[i] == (b'G'){
            bytes[i] = b'C';
        }
        else if bytes[i] == (b'T'){
            bytes[i] = b'A';
        }
        else { 
            //bytes[i] = b'X';
        }
    }
    println!("{}", std::str::from_utf8(bytes.as_slice()).expect(""));
}
