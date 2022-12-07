mod histogram;

use std::env;
use std::fs;
use std::io;
use std::io::BufRead;
use std::io::Write;

use histogram::print_histogram;

fn main(){
    let args: Vec<String> = env::args().collect();
    let filename = args[1].clone();

    // Histogram parameters
    let min: i64 = args[2].parse().unwrap();
    let max: i64 = args[3].parse().unwrap();
    let nbins: i64 = args[4].parse().unwrap();

    let file = fs::File::open(filename).unwrap();
    let lines = io::BufReader::new(file).lines();

    let mut starts = Vec::<i64>::new();
    for line in lines{
        let line = line.unwrap();
        let tokens: Vec<&str> = line.split('\t').collect();
        let length: u64 = tokens[2].parse().unwrap();
        starts.push(length as i64);
    }

    print_histogram(starts.iter().map(|x| *x), min, max, nbins);
}