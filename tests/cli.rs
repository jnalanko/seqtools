
use std::process::{Command, Stdio}; // Run programs
use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use tempfile;
use std::str;
use std::io::Write;

#[test]
fn stats() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("seqtools")?;

    cmd.arg("stats").arg("tests/data/reads.fastq.gz");
    cmd.assert()
        .stdout(predicate::str::contains("Number of nucleotides: 474"))
        .stdout(predicate::str::contains("Number of sequences: 10"))
        .stdout(predicate::str::contains("Maximum sequence length: 54"))
        .stdout(predicate::str::contains("Minimum sequence length: 38"))
        .stdout(predicate::str::contains("Average sequence length: 47.4"));

    Ok(())
}

#[test]
fn convert() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("seqtools")?;

    // fasta -> fasta.gz
    let first = cmd.arg("convert").arg("tests/data/reads.fna").arg("--fasta-out").arg("--gzip-out").stdout(Stdio::piped()).spawn()?;
    let first_stdout_result = first.wait_with_output()?.stdout;

    // fasta.gz -> fastq
    cmd = Command::cargo_bin("seqtools")?;
    let mut second = cmd.arg("convert").arg("--fasta-in").arg("--gzip-in").arg("--fastq-out").stdin(Stdio::piped()).stdout(Stdio::piped()).spawn()?;
    let second_stdin = second.stdin.as_mut().unwrap();
    second_stdin.write_all(&first_stdout_result)?;
    drop(second_stdin); // Flush
    let second_stdout_result = second.wait_with_output()?.stdout;

    // fastq -> fastq.gz
    cmd = Command::cargo_bin("seqtools")?;
    let mut third = cmd.arg("convert").arg("--fastq-in").arg("--fastq-out").arg("--gzip-out").stdin(Stdio::piped()).stdout(Stdio::piped()).spawn()?;
    let third_stdin = third.stdin.as_mut().unwrap();
    third_stdin.write_all(&second_stdout_result)?;
    drop(third_stdin); // Flush
    let third_stdout_result = third.wait_with_output()?.stdout;

    // fastq.gz -> fasta
    cmd = Command::cargo_bin("seqtools")?;
    let mut fourth = cmd.arg("convert").arg("--fastq-in").arg("--gzip-in").arg("--fasta-out").stdin(Stdio::piped()).stdout(Stdio::piped()).spawn()?;
    let fourth_stdin = fourth.stdin.as_mut().unwrap();
    fourth_stdin.write_all(&third_stdout_result)?;
    drop(fourth_stdin); // Flush
    let fourth_stdout_result = fourth.wait_with_output()?.stdout;

    // See if we have the same fasta data as we started with
    let S = str::from_utf8(&fourth_stdout_result)?;
    let answer = "\
>SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1
TTGGACCGGCGCAAGACGGACCAGNGCGAAAGCATTTGCCAAGAANNNN
>SRR403017.2 HWUSI-EAS108E_0007:3:1:10327:976/1
CAACTTTCTATCTGGCATTCCCTGNGGAGGAAATAGAATCGCNNNN
>SRR403017.3 HWUSI-EAS108E_0007:3:1:13569:972/1
GATCGGAAGAGCACACGTCTGAACNCCAGTCACTTAGGCATCTCGNNNN
>SRR403017.4 HWUSI-EAS108E_0007:3:1:16652:968/1
TCTGTTATGTCGGATTACTGCTGTNAGTCAGTGNNNNN
>SRR403017.5 HWUSI-EAS108E_0007:3:1:6652:978/1
CTGTGTTTTTATTTAGTGGGTAGCNGAAGTTGTTCAGAAGAGCAGNNNN
>SRR403017.6 HWUSI-EAS108E_0007:3:1:8133:984/1
TGGAGCCCTTAAAAACGTACAGGCNCCTGGGGTCATTGGGGGTGAGGTNNNN
>SRR403017.7 HWUSI-EAS108E_0007:3:1:8545:977/1
TTAGATGTCCGGGGCTGCACGTGCNCTATGACTGGCTCAGCGNNNN
>SRR403017.8 HWUSI-EAS108E_0007:3:1:8825:983/1
ATAGAGAAGGGGGACAATGAGCCTGGATCTTTGCCTTGNNNN
>SRR403017.9 HWUSI-EAS108E_0007:3:1:11308:980/1
AAATTTATTTCTCATAGTTCTGGANTCTAGGAAGTTCAAGATCAGNNNN
>SRR403017.10 HWUSI-EAS108E_0007:3:1:14685:981/1
CCCATTCTTGGAGATACCAGCAAAAATTCNAATTCACCAACACCAGCAGCNNNN
";

    assert_eq!(S, answer);
    Ok(())
}



#[test]
fn trim() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("seqtools")?;

    cmd.arg("trim").arg("tests/data/reads.fastq.gz")
        .arg("--from-start").arg("20").arg("--from-end").arg("26").arg("--fastq-out");
    // Corner case: 20 + 26 = 46. There is a sequence of length exactly 46. That should be deleted.

    let answer = 
    "\
@SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1
CCA
+
567
@SRR403017.3 HWUSI-EAS108E_0007:3:1:13569:972/1
GAA
+
III
@SRR403017.5 HWUSI-EAS108E_0007:3:1:6652:978/1
TAG
+
III
@SRR403017.6 HWUSI-EAS108E_0007:3:1:8133:984/1
AGGCNC
+
IIIIII
@SRR403017.9 HWUSI-EAS108E_0007:3:1:11308:980/1
TGG
+
III
@SRR403017.10 HWUSI-EAS108E_0007:3:1:14685:981/1
CAAAAATT
+
IIIIIIII";

    let deleted_seqs_message = "Deleted 4 sequences as too short to trim";

    cmd.assert()
        .stdout(predicate::str::contains(answer))
        .stderr(predicate::str::contains(deleted_seqs_message));

    Ok(())
}