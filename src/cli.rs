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
            Arg::new("fasta")
                .short('a')
                .long("fasta")
                .action(ArgAction::SetTrue)
                .help("Parse in fasta format")
                .global(true),
        )
        .arg(
            Arg::new("fastq")
                .short('q')
                .long("fastq")
                .action(ArgAction::SetTrue)
                .help("Parse in fastq format")
                .global(true),
        )
        .arg(
            Arg::new("gzip")
                .short('g')
                .long("gzip")
                .action(ArgAction::SetTrue)
                .help("Read gzipped input")
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
                ),
        )
        .subcommand(Command::new("stats").about("Print stats about the input."))
}
