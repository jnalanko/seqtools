
use std::process::{Command, Stdio}; // Run programs
use assert_cmd::prelude::*; 
// Add methods on commands
use predicates::prelude::*; // Used for writing assertions

use std::str;
use std::io::{Write, BufWriter};

#[test]
fn stats() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("seqtools")?;

    cmd.arg("stats").arg("tests/data/reads.fastq.gz");
    cmd.assert()
        .stdout(predicate::str::contains("Number of nucleotides: 474"))
        .stdout(predicate::str::contains("Number of sequences: 10"))
        .stdout(predicate::str::contains("Maximum sequence length: 54"))
        .stdout(predicate::str::contains("Minimum sequence length: 38"))
        .stdout(predicate::str::contains("Average sequence length: 47.4"))
        .stdout(predicate::str::contains("Maximum quality value: 93"))
        .stdout(predicate::str::contains("Minimum quality value: 0"))
        .stdout(predicate::str::contains("Average quality value: 41.4008"));
    Ok(())
}

#[test]
fn length_histogram() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("seqtools")?;

    cmd.arg("length-histogram").arg("tests/data/reads.fna").arg("--min").arg("10").arg("--max").arg("100").arg("--nbins").arg("15");

    let answer = "\
10	
16	
22	
28	
34	######
40	######
46	########################################
52	#############
58	
64	
70	
76	
82	
88	
94	
";
    cmd.assert()
        .stdout(predicate::str::contains(answer));

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
    let mut second = cmd.arg("convert").arg("--gzip-in").arg("--fastq-out").stdin(Stdio::piped()).stdout(Stdio::piped()).spawn()?;
    let second_stdin = second.stdin.as_mut().unwrap();
    second_stdin.write_all(&first_stdout_result)?;
    let second_stdout_result = second.wait_with_output()?.stdout;

    // fastq -> fastq.gz
    cmd = Command::cargo_bin("seqtools")?;
    let mut third = cmd.arg("convert").arg("--fastq-out").arg("--gzip-out").stdin(Stdio::piped()).stdout(Stdio::piped()).spawn()?;
    let third_stdin = third.stdin.as_mut().unwrap();
    third_stdin.write_all(&second_stdout_result)?;
    let third_stdout_result = third.wait_with_output()?.stdout;

    // fastq.gz -> fasta
    cmd = Command::cargo_bin("seqtools")?;
    let mut fourth = cmd.arg("convert").arg("--gzip-in").arg("--fasta-out").stdin(Stdio::piped()).stdout(Stdio::piped()).spawn()?;
    let fourth_stdin = fourth.stdin.as_mut().unwrap();
    fourth_stdin.write_all(&third_stdout_result)?;
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
fn subsample_frac() -> Result<(), Box<dyn std::error::Error>>{
    // Test fraction subsampling
    let mut cmd = Command::cargo_bin("seqtools")?;
    let child = cmd.arg("subsample").arg("tests/data/reads.fna").arg("--fraction").arg("0.55").arg("--fasta-out").stdout(Stdio::piped()).spawn()?;
    let child_out = child.wait_with_output()?.stdout;

    let n_lines = child_out.iter().filter(|&&c| c == b'\n').count();
    assert_eq!(n_lines, 5*2); // 5 sequences and headers

    Ok(())
}

#[test]
fn subsample_howmany() -> Result<(), Box<dyn std::error::Error>>{
    // Test howmany subsampling
    let mut cmd = Command::cargo_bin("seqtools")?;
    let child = cmd.arg("subsample").arg("tests/data/reads.fna").arg("--howmany").arg("4").arg("--fasta-out").stdout(Stdio::piped()).spawn()?;
    let child_out = child.wait_with_output()?.stdout;

    let n_lines = child_out.iter().filter(|&&c| c == b'\n').count();
    assert_eq!(n_lines, 4*2); // 4 sequences and headers

    Ok(())
}

#[test]
fn subsample_toomany() -> Result<(), Box<dyn std::error::Error>>{
    // Test subsampling more than the number of sequences in the file
    let mut cmd = Command::cargo_bin("seqtools")?;
    let child = cmd.arg("subsample").arg("tests/data/reads.fna").arg("--howmany").arg("11").arg("--fasta-out").stdout(Stdio::piped()).spawn()?;
    let child_out = child.wait_with_output()?.stdout;

    let n_lines = child_out.iter().filter(|&&c| c == b'\n').count();
    assert_eq!(n_lines, 10*2); // 10 sequences and headers

    Ok(())
}


#[test]
fn subsample() -> Result<(), Box<dyn std::error::Error>>{

    subsample_frac()?;
    subsample_howmany()?;
    subsample_toomany()?;

    Ok(())
}

#[test]
fn remove_duplicates() -> Result<(), Box<dyn std::error::Error>>{
    let buf = Vec::<u8>::new();
    let mut bufwriter = BufWriter::<Vec::<u8>>::new(buf);
    let mut seqwriter = jseqio::writer::FastXWriter::new(&mut bufwriter, jseqio::FileType::FASTA);
    seqwriter.write(&jseqio::record::OwnedRecord{
        head: b"".to_vec(),
        seq: b"ACGT".to_vec(),
        qual: None,
    });
    seqwriter.write(&jseqio::record::OwnedRecord{
        head: b"".to_vec(),
        seq: b"GGG".to_vec(),
        qual: None,
    });
    seqwriter.write(&jseqio::record::OwnedRecord{
        head: b"".to_vec(),
        seq: b"ACGT".to_vec(),
        qual: None,
    });
    seqwriter.write(&jseqio::record::OwnedRecord{
        head: b"".to_vec(),
        seq: b"ACGT".to_vec(),
        qual: None,
    });

    seqwriter.flush();
    drop(seqwriter);
    let buf = bufwriter.into_inner().unwrap(); // Take back ownership

    let mut cmd = Command::cargo_bin("seqtools")?;
    let mut child = cmd.arg("remove-duplicates").arg("--fasta-out").stdout(Stdio::piped()).stdin(Stdio::piped()).spawn()?;
    child.stdin.take().unwrap().write_all(buf.as_slice()).unwrap();
    // TODO: Do we need to write eof?
    let child_out = child.wait_with_output()?.stdout;

    assert_eq!(child_out, b">\nACGT\n>\nGGG\n");

    Ok(())
}

#[test]
fn trim() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("seqtools")?;

    cmd.arg("trim").arg("tests/data/reads.fastq.gz")
        .arg("--from-start").arg("20").arg("--from-end").arg("26").arg("--fastq-out").arg("--min-final-length").arg("4");

    let answer = 
    "\
@SRR403017.6 HWUSI-EAS108E_0007:3:1:8133:984/1
AGGCNC
+
IIIIII
@SRR403017.10 HWUSI-EAS108E_0007:3:1:14685:981/1
CAAAAATT
+
IIIIIIII";

    let deleted_seqs_message = "Deleted 8 sequences whose final length would have been below the minimum length 4";

    cmd.assert()
        .stdout(predicate::str::contains(answer))
        .stderr(predicate::str::contains(deleted_seqs_message));

    Ok(())
}