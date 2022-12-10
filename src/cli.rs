use clap::{Arg, ArgAction, Command};

pub fn build_cli() -> Command {
    Command::new("seq_tools")
        .version("0.1.0")
        .author("Jarno N. Alanko <alanko.jarno@gmail.com>")
        .about("Fasta/fastq parsing")
        .arg_required_else_help(true)
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .help("Input filename")
                .global(true),
        )
        .arg(
            Arg::new("fasta-in")
                .long("fasta-in")
                .action(ArgAction::SetTrue)
                .help("Parse in fasta format")
                .global(true),
        )
        .arg(
            Arg::new("fastq-in")
                .long("fastq-in")
                .action(ArgAction::SetTrue)
                .help("Parse in fastq format")
                .global(true),
        )
        .arg(
            Arg::new("gzip-in")
                .long("gzip-in")
                .action(ArgAction::SetTrue)
                .help("Read gzipped input")
                .global(true),
        )
        .arg(
            Arg::new("fasta-out")
                .long("fasta-out")
                .action(ArgAction::SetTrue)
                .help("Write in fasta format")
                .global(true),
        )
        .arg(
            Arg::new("fastq-out")
                .long("fastq-out")
                .action(ArgAction::SetTrue)
                .help("Write in fastq format")
                .global(true),
        )
        .arg(
            Arg::new("gzip-out")
                .long("gzip-out")
                .action(ArgAction::SetTrue)
                .help("Write as gzipped")
                .global(true),
        )
        .subcommand(
            Command::new("length-histogram")
                .about("Print the length histogram of the sequences.")
                .arg_required_else_help(true)
                .arg(Arg::new("min").default_value("0").help("Minimum value"))
                .arg(Arg::new("max").default_value("1000").help("Maximum value"))
                .arg(Arg::new("nbins").default_value("20").help("Number of bins")),
        ).subcommand(
            Command::new("subsample")
                .about("Subsample a random fraction of sequences.")
                .arg_required_else_help(true)
                .arg(Arg::new("fraction")
                    .short('f')
                    .long("fraction")
                    .required(true)
                ).arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Output filename")
                        .global(true),
                ),
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
                ).arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Output filename")
                        .global(true),
                ),
        ).subcommand(
            Command::new("convert")
                .about("Convert the input file format into the output file format.")
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Output filename")
                        .global(true),
                ),
        )
        .subcommand(Command::new("stats").about("Print stats about the input."))
}
