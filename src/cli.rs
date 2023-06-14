use clap::{Arg, ArgAction, Command};

pub fn build_cli() -> Command {
    
    let stdout_fasta = 
        Arg::new("fasta-out")
            .long("fasta-out")
            .action(ArgAction::SetTrue)
            .help("Write to stdout in fasta format")
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

    Command::new("seq_tools")
        .version("0.1.0")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .arg_required_else_help(true)
        .arg(
            Arg::new("input")
                //.short('i')
                //.long("input")
                .help("Read from a file on disk (format is deduced from the file extension)")
                .global(true),
        )
        .arg(
            Arg::new("fasta-in")
                .long("fasta-in")
                .action(ArgAction::SetTrue)
                .help("Parse from stdin in fasta format")
                .global(true),
        )
        .arg(
            Arg::new("fastq-in")
                .long("fastq-in")
                .action(ArgAction::SetTrue)
                .help("Parse from stdin in fastq format")
                .global(true),
        )
        .arg(
            Arg::new("gzip-in")
                .long("gzip-in")
                .action(ArgAction::SetTrue)
                .help("Parse from stdin as gzipped")
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
                .about("Subsample a random fraction of sequences.")
                .arg_required_else_help(true)
                .arg(Arg::new("fraction")
                    .short('f')
                    .long("fraction")
                    .required(true)
                    .conflicts_with("howmany")
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
            Command::new("extract-read")
                .about("Prints the read with rank i (zero-based)")
                .arg_required_else_help(true)                
                .arg(Arg::new("rank")
                    .short('r')
                    .long("rank")
                    .required(true)
                )                
        )
        .subcommand(Command::new("stats").about("Print stats about the input."))
}
