extern crate flate2;

use seq_tools::*;
use trim_adapters::TrimMode;

mod cli;

fn read_lines(filename: &str) -> Vec<String>{
    let reader = std::io::BufReader::new(std::fs::File::open(filename).unwrap());
    std::io::BufRead::lines(reader).map(|s| s.unwrap().to_owned()).collect()
}


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
        Some(("print-lengths", sub_matches)) => { 
            let mut reader = get_reader(&matches).unwrap();
            print_lengths(&mut reader);
        }
        Some(("stats", _)) => { 
            let mut reader = get_reader(&matches).unwrap();
            print_stats(&mut reader);
        }
        Some(("extract-reads", sub_matches)) => { 
            let reader = get_reader(&matches).unwrap();
            if let Some(ranks) = sub_matches.get_many::<String>("rank"){
                let list: Vec<usize> = ranks.map(|s| s.parse::<usize>().unwrap()).collect();
                extract_reads_by_ranks(reader, &list);
            } else if let Some(ranks_listfilename) = sub_matches.get_one::<String>("ranks-listfile"){
                let ranks_reader = std::io::BufReader::new(std::fs::File::open(ranks_listfilename).unwrap());
                let list = std::io::BufRead::lines(ranks_reader).map(|s| s.unwrap().parse::<usize>().unwrap()).collect();
                extract_reads_by_ranks(reader, &list);
            } else if let Some(names) = sub_matches.get_many::<String>("name"){
                let list: Vec<String> = names.map(|s| s.to_owned()).collect();
                extract_reads_by_names(reader, &list);
            } else if let Some(names_listfilename) = sub_matches.get_one::<String>("names-listfile"){
                let list = read_lines(names_listfilename);
                extract_reads_by_names(reader, &list);
            }
        }
        Some(("subsample", sub_matches)) => { // TODO: Untested
            if matches.get_one::<String>("input").is_none() {
                panic!("Can not subsample from stdin because we need to pass over the data twice.");
            }

            let seed = sub_matches.get_one::<String>("seed").map(|s| s.parse::<u64>().unwrap());
            let paired_interleaved = sub_matches.get_flag("paired-interleaved");

            if let Some(f) = sub_matches.get_one::<String>("fraction"){
                let frac = f.parse::<f64>().unwrap();
                // Get two readers for two passes over the data
                let input1 = get_reader(&matches).unwrap();
                let input2 = get_reader(&matches).unwrap();
                let mut output = get_writer(sub_matches);
                random_subsample(input1,  input2, &mut output, frac, seed, paired_interleaved);
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
                random_subsample_howmany(input, &mut output, total_seqs as usize, howmany as usize, seed, paired_interleaved);
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
        Some(("reverse-complement", sub_matches)) => { 
            let mut reader = get_reader(&matches).unwrap();
            let mut writer = get_writer(sub_matches);
            reverse_complement(&mut reader, &mut writer);
        }
        Some(("concat", sub_matches)) => { 
            let header = match sub_matches.get_one::<String>("header"){
                Some(header) => header.as_bytes().to_owned(),
                None => "".as_bytes().to_owned(),
            };
            let mut reader = get_reader(&matches).unwrap();
            let mut writer = get_writer(sub_matches);
            concatenate(&mut reader, &mut writer, &header);
        }
        Some(("trim", sub_matches)) => { 
            let mut reader = get_reader(&matches).unwrap();
            let mut writer = get_writer(sub_matches);
            let from_start: usize = sub_matches.get_one::<String>("from-start").unwrap().parse().unwrap();
            let from_end: usize = sub_matches.get_one::<String>("from-end").unwrap().parse().unwrap();
            let min_final_length: usize = sub_matches.get_one::<String>("min-final-length").unwrap().parse().unwrap();
            trim(&mut reader, &mut writer, from_start, from_end, min_final_length);
        }
        Some(("trim-adapters", sub_matches)) => { 
            let mut reader = get_reader(&matches).unwrap();
            let mut writer = get_writer(sub_matches);
            let adapter_file = sub_matches.get_one::<std::path::PathBuf>("adapters").unwrap();
            let min_final_length: usize = sub_matches.get_one::<String>("min-final-length").unwrap().parse().unwrap();
            let max_trim_length: usize = sub_matches.get_one::<String>("max-trim-length").unwrap().parse().unwrap();
            let identity_threshold = *sub_matches.get_one::<f64>("identity-threshold").unwrap();

            let adapter_lines: Vec<String> = read_lines(adapter_file.to_str().unwrap());
            let adapters: Vec<(Vec<u8>, TrimMode)> = adapter_lines.into_iter().map(|line| {
                let mut parts = line.split_whitespace();
                assert_eq!(parts.clone().count(), 2);
                let mode_string = parts.next().unwrap();
                let mode = match mode_string {
                    "upto" => TrimMode::Upto,
                    "from" => TrimMode::From,
                    _ => panic!("Invalid trim mode {}: should be 'upto' or 'from'", mode_string),
                };
                let adapter = parts.next().unwrap().as_bytes().to_vec();
                assert!(adapter.is_ascii());
                (adapter, mode)
            }).collect();

            // Print loaded adapters and their trim modes 
            for (adapter, trim_mode) in adapters.iter(){
                eprintln!("Loaded adapter: {}, Trim mode: {:?}", String::from_utf8_lossy(adapter), trim_mode);
                for &c in adapter{
                    assert!(c == b'A' || c == b'C' || c == b'G' || c == b'T' || c == b'N');
                }
            }

            seq_tools::trim_adapters::trim_adapters(&mut reader, &mut writer, adapters, max_trim_length, min_final_length, identity_threshold);
        }
        _ => {}
    };
}
