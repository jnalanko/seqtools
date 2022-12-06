use std::env;
use std::fs;
use std::io;
use std::io::BufRead;
use std::io::Write;

fn print_histogram(numbers: &Vec<u64>, min: i64, max: i64, n_bins: i64){

    let mut counters: Vec<i64> = vec![0; n_bins as usize];
    let bin_width = (max-min+1) / n_bins;
    for len in numbers{
        let mut bin = (*len as i64 - min as i64) / bin_width;

        // Clamp to [0, n_bins-1]
        bin = std::cmp::max(0, bin);
        bin = std::cmp::min(n_bins-1, bin);

        counters[bin as usize] += 1;
    }

    let max_counter: i64 = *counters.iter().max().unwrap();
    let n_columns: i64 = 40;

    for (i, c) in counters.iter().enumerate(){
        let n_chars = ((*c as f64 / max_counter as f64) * n_columns as f64) as i64;
        print!("{}\t", (min + (i as i64)*bin_width) as usize);
        std::io::stdout().write_all(vec![b'#'; n_chars as usize].as_slice()).ok();
        println!();
    }
}

fn main(){
    let args: Vec<String> = env::args().collect();
    let filename = args[1].clone();

    // Histogram parameters
    let min: u64 = args[2].parse().unwrap();
    let max: u64 = args[3].parse().unwrap();
    let nbins: u64 = args[4].parse().unwrap();

    let file = fs::File::open(filename).unwrap();
    let lines = io::BufReader::new(file).lines();

    let mut starts = Vec::<u64>::new();
    for line in lines{
        let line = line.unwrap();
        let tokens: Vec<&str> = line.split('\t').collect();
        let length: u64 = tokens[2].parse().unwrap();
        starts.push(length);
    }

    print_histogram(&starts, min as i64, max as i64, nbins as i64);
}