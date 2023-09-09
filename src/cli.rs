use clap::{Arg, ArgAction, Command};

pub fn build_cli() -> Command {
    
    let stdout_fasta = 
        Arg::new("fasta-out")
            .long("fasta-out")
            .action(ArgAction::SetTrue)
            .help("Write to stdout in fasta format")
            .conflicts_with("fastq-out")
            .global(false);
    let stdout_fastq =
        Arg::new("fastq-out")
            .long("fastq-out")
            .action(ArgAction::SetTrue)
            .help("Write to stdout in fastq format")
            .global(false);
    let stdout_gzip =
        Arg::new("gzip-out")
            .long("gzip-out")
            .action(ArgAction::SetTrue)
            .help("Write to stdout as gzipped")
            .global(false);
    let output_file =
        Arg::new("output")
        .short('o')
        .long("output")
        .help("Output filename")
        .global(false);

    Command::new("seqtools")
        .version("0.1.0")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .arg_required_else_help(true)
        .arg(
            Arg::new("input")
                //.short('i')
                //.long("input")
                .help("Read from a file on disk. Format can be fasta or fastq, gzipped or not. The format is detected automatically.)")
                .global(true),
        )
        .subcommand(
            Command::new("length-histogram")
                .about("Print the length histogram of the sequences.")
                .arg_required_else_help(true)
                .arg(Arg::new("min").long("min").default_value("0").help("Minimum value"))
                .arg(Arg::new("max").long("max").default_value("1000").help("Maximum value"))
                .arg(Arg::new("nbins").long("nbins").default_value("20").help("Number of bins")),
        ).subcommand(
            Command::new("subsample")
                .about("Randomly subsample a given fraction of sequences.")
                .arg_required_else_help(true)
                .arg(Arg::new("fraction")
                    .short('f')
                    .long("fraction")
                    .required(true)
                    .conflicts_with("howmany")
                )
                .arg(Arg::new("seed")
                    .help("The seed for the random number generator. If not given, the seed is generated from the current time.")
                    .short('s')
                    .long("seed")
                    .required(false)
                )
                .arg(Arg::new("howmany")
                    .short('n')
                    .long("howmany")
                    .required(true)
                    .conflicts_with("fraction")
                )
                .arg(&output_file)
                .arg(&stdout_fasta)
                .arg(&stdout_fastq)
                .arg(&stdout_gzip),
        ).subcommand(
            Command::new("remove-duplicates")
                .about("Removes reads that have exactly the same nucleotides (headers do not need to match).")
                .arg_required_else_help(true)
                .arg(&output_file)
                .arg(&stdout_fasta)
                .arg(&stdout_fastq)
                .arg(&stdout_gzip),
            ).subcommand(
            Command::new("trim")
                .about("Trim starts and ends of sequences.")
                .arg_required_else_help(true)
                .arg(Arg::new("from-start")
                    .long("from-start")
                    .required(true)
                ).arg(Arg::new("from-end")
                    .long("from-end")
                    .required(true)
                ).arg(Arg::new("min-final-length")
                    .long("min-final-length")
                    .default_value("1")
                    .required(true)
                ).arg(&output_file)
                .arg(&stdout_fasta)
                .arg(&stdout_fastq)
                .arg(&stdout_gzip),
        ).subcommand(
            Command::new("convert")
                .about("Convert the input file format into the output file format.")
                .arg_required_else_help(true)                
                .arg(&output_file)
                .arg(&stdout_fasta)
                .arg(&stdout_fastq)
                .arg(&stdout_gzip),
        ).subcommand(
            Command::new("extract-reads")
                .about("Extract read by rank or sequence name.")
                .arg_required_else_help(true)                
                .arg(Arg::new("rank")
                    .short('r')
                    .long("rank")
                    .help("Zero-based ranks. Can be given multiple times (e.g. -r 4 -r 6).")
                    .conflicts_with("name").conflicts_with("names-listfile")
                    .action(ArgAction::Append) // Can have multiple
                )
                .arg(Arg::new("ranks-listfile")
                    .long("ranks-listfile")
                    .help("Zero-based ranks. One number per line.")
                    .conflicts_with("name").conflicts_with("names-listfile")
                )
                .arg(Arg::new("name")
                    .short('a')
                    .long("name")
                    .help("Can be given multiple times. Name is the part of the header that comes before the first space character, without the leading '>' or '@'.")
                    .action(ArgAction::Append) // Can have multiple
                )
                    .arg(Arg::new("names-listfile")
                    .long("names-listfile")
                    .help("One sequence name per line. Name is the part of the header that comes before the first space character, without the leading '>' or '@'.")
                )
        )
        .subcommand(Command::new("stats").about("Print stats about the input."))
}
