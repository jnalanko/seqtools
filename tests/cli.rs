
use std::process::Command; // Run programs
use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions

#[test]
fn test_stats() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("seq_tools")?;

    cmd.arg("stats").arg("-i").arg("tests/data/reads.fastq.gz");
    cmd.assert()
        .stdout(predicate::str::contains("Number of nucleotides: 474"))
        .stdout(predicate::str::contains("Number of sequences: 10"))
        .stdout(predicate::str::contains("Maximum sequence length: 54"))
        .stdout(predicate::str::contains("Minimum sequence length: 38"))
        .stdout(predicate::str::contains("Average sequence length: 47.4"));

    Ok(())
}

fn test_trim() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("seq_tools")?;

    cmd.arg("trim").arg("-i").arg("tests/data/reads.fastq.gz")
        .arg("from-start").arg("20").arg("from-end").arg("22").arg("--fastq-out");

    let answer = "@SRR403017.1 HWUSI-EAS108E_0007:3:1:3797:973/1
    CCAGNGC
    +
    56789:;
    @SRR403017.2 HWUSI-EAS108E_0007:3:1:10327:976/1
    CCTG
    +
    fghi
    @SRR403017.3 HWUSI-EAS108E_0007:3:1:13569:972/1
    GAACNCC
    +
    IIIIIII
    @SRR403017.5 HWUSI-EAS108E_0007:3:1:6652:978/1
    TAGCNGA
    +
    IIIIIII
    @SRR403017.6 HWUSI-EAS108E_0007:3:1:8133:984/1
    AGGCNCCTGG
    +
    IIIIIIIIII
    @SRR403017.7 HWUSI-EAS108E_0007:3:1:8545:977/1
    GTGC
    +
    IIII
    @SRR403017.9 HWUSI-EAS108E_0007:3:1:11308:980/1
    TGGANTC
    +
    IIIIIII
    @SRR403017.10 HWUSI-EAS108E_0007:3:1:14685:981/1
    CAAAAATTCNAA
    +
    IIIIIIIIIIII
    ";

    let deleted_seqs_message = "Deleted 2 sequences as too short to trim";

    cmd.assert()
        .stdout(predicate::str::contains(answer))
        .stdout(predicate::str::contains(deleted_seqs_message));

    Ok(())
}