extern crate flate2;

use my_seqio::reader::DynamicFastXReader;
use my_seqio::writer::DynamicFastXWriter;

use std::io::{Write,BufWriter};
use rand::Rng;
use std::cmp::{max,min};
use seq_tools::*;

mod cli;

fn main() {

    let matches = cli::build_cli().get_matches();

    match matches.subcommand() {
        Some(("length-histogram", sub_matches)) => { 
            let mut reader = get_reader(&matches);
            let min: i64 = sub_matches.get_one::<String>("min").unwrap().parse::<i64>().unwrap();
            let max: i64 = sub_matches.get_one::<String>("max").unwrap().parse::<i64>().unwrap();
            let nbins: i64 = sub_matches.get_one::<String>("nbins").unwrap().parse::<i64>().unwrap();
            print_length_histogram(&mut reader, min as i64, max as i64, nbins as i64);
        }
        Some(("stats", _)) => { 
            let mut reader = get_reader(&matches);
            print_stats(&mut reader);
        }
        Some(("extract-read", sub_matches)) => { 
            let mut reader = get_reader(&matches);
            let rank = sub_matches.get_one::<String>("rank").unwrap().parse::<usize>().unwrap();
            extract_read(&mut reader, rank);
        }
        Some(("subsample", sub_matches)) => { // TODO: Untested
            if matches.get_one::<String>("input") == None {
                panic!("Can not subsample from stdin because we need to pass over the data twice.");
            }

            if let Some(f) = sub_matches.get_one::<String>("fraction"){
                let frac = f.parse::<f64>().unwrap();
                // Get two readers for two passes over the data
                let input1 = get_reader(&matches);
                let input2 = get_reader(&matches);
                let mut output = get_writer(&sub_matches);
                random_subsample(input1,  input2, &mut output, frac);
            }

            if let Some(f) = sub_matches.get_one::<String>("howmany"){
                let mut howmany = f.parse::<u64>().unwrap();

                // Count the number of sequences in the file
                eprintln!("Counting sequences...");
                let total_seqs = count_sequences(get_reader(&matches));
                eprintln!("{} sequences found", total_seqs);
                eprintln!("Subsampling {} sequences...", howmany);
                if howmany > total_seqs{
                    eprintln!("Warning: Trying to sample more sequences than what the file has -> sampling all.");
                    howmany = total_seqs;
                }

                // Do the subsampling
                let input = get_reader(&matches);
                let mut output = get_writer(&sub_matches);
                random_subsample_howmany(input, &mut output, total_seqs as usize, howmany as usize);
            }
        }
        Some(("convert", sub_matches)) => { 
            let mut reader = get_reader(&matches);
            let mut writer = get_writer(&sub_matches);
            convert(&mut reader, &mut writer);
        }
        Some(("trim", sub_matches)) => { 
            let mut reader = get_reader(&matches);
            let mut writer = get_writer(&sub_matches);
            let from_start: usize = sub_matches.get_one::<String>("from-start").unwrap().parse().unwrap();
            let from_end: usize = sub_matches.get_one::<String>("from-end").unwrap().parse().unwrap();
            trim(&mut reader, &mut writer, from_start, from_end);
        }
        _ => {}
    };
}
