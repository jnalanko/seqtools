use std::io::Write;

// Takes an iterator that produces i64 values, and prints the histogram
// of those values to stdout.
fn print_histogram(value_iterator: impl Iterator<Item = i64>, min: i64, max: i64, n_bins: i64){
    let mut counters: Vec<i64> = vec![0; n_bins as usize];
    let bin_width = (max-min+1) / n_bins;

    for x in value_iterator{
        let mut bin = (x - min) / bin_width;

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
