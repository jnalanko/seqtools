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
    let bytes = args[1].as_bytes();
    let n = bytes.len();
    for i in 0..bytes.len(){
        if bytes[n-1-i] == ('A' as u8){print!("{}",'T')};
        if bytes[n-1-i] == ('C' as u8){print!("{}",'G')};
        if bytes[n-1-i] == ('G' as u8){print!("{}",'C')};
        if bytes[n-1-i] == ('T' as u8){print!("{}",'A')};
    }
    println!();
}
