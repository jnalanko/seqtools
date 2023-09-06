extern crate flate2;

use seq_tools::*;

mod cli;

fn main() {

    let matches = cli::build_cli().get_matches();

    match matches.subcommand() {
        Some(("length-histogram", sub_matches)) => { 
            let mut reader = get_reader(&matches).unwrap();
            let min: i64 = sub_matches.get_one::<String>("min").unwrap().parse::<i64>().unwrap();
            let max: i64 = sub_matches.get_one::<String>("max").unwrap().parse::<i64>().unwrap();
            let nbins: i64 = sub_matches.get_one::<String>("nbins").unwrap().parse::<i64>().unwrap();
            print_length_histogram(&mut reader, min, max, nbins);
        }
        Some(("stats", _)) => { 
            let mut reader = get_reader(&matches).unwrap();
            print_stats(&mut reader);
        }
        Some(("extract-reads", sub_matches)) => { 
            let mut reader = get_reader(&matches).unwrap();
            let ranks: Vec<usize> = sub_matches.get_many::<String>("rank").unwrap().map(|s| s.parse::<usize>().unwrap()).collect();
            extract_reads(&mut reader, &ranks); //todo: header prefix
        }
        Some(("subsample", sub_matches)) => { // TODO: Untested
            if matches.get_one::<String>("input").is_none() {
                panic!("Can not subsample from stdin because we need to pass over the data twice.");
            }

            let seed = if let Some(f) = sub_matches.get_one::<String>("seed"){
                Some(f.parse::<u64>().unwrap())
            } else{
                None
            };

            if let Some(f) = sub_matches.get_one::<String>("fraction"){
                let frac = f.parse::<f64>().unwrap();
                // Get two readers for two passes over the data
                let input1 = get_reader(&matches).unwrap();
                let input2 = get_reader(&matches).unwrap();
                let mut output = get_writer(sub_matches);
                random_subsample(input1,  input2, &mut output, frac, seed);
            }

            if let Some(f) = sub_matches.get_one::<String>("howmany"){
                let mut howmany = f.parse::<u64>().unwrap();

                // Count the number of sequences in the file
                eprintln!("Counting sequences...");
                let total_seqs = count_sequences(get_reader(&matches).unwrap());
                eprintln!("{} sequences found", total_seqs);
                eprintln!("Subsampling {} sequences...", howmany);
                if howmany > total_seqs{
                    eprintln!("Warning: Trying to sample more sequences than what the file has -> sampling all.");
                    howmany = total_seqs;
                }

                // Do the subsampling
                let input = get_reader(&matches).unwrap();
                let mut output = get_writer(sub_matches);
                random_subsample_howmany(input, &mut output, total_seqs as usize, howmany as usize, seed);
            }
        }
        Some(("remove-duplicates", sub_matches)) => { // TODO: Untested

            let mut input = get_reader(&matches).unwrap();
            let mut output = get_writer(sub_matches);
            remove_duplicates(&mut input, &mut output);
        }
        Some(("convert", sub_matches)) => { 
            let mut reader = get_reader(&matches).unwrap();
            let mut writer = get_writer(sub_matches);
            convert(&mut reader, &mut writer);
        }
        Some(("trim", sub_matches)) => { 
            let mut reader = get_reader(&matches).unwrap();
            let mut writer = get_writer(sub_matches);
            let from_start: usize = sub_matches.get_one::<String>("from-start").unwrap().parse().unwrap();
            let from_end: usize = sub_matches.get_one::<String>("from-end").unwrap().parse().unwrap();
            let min_final_length: usize = sub_matches.get_one::<String>("min-final-length").unwrap().parse().unwrap();
            trim(&mut reader, &mut writer, from_start, from_end, min_final_length);
        }
        _ => {}
    };
}
