
use std::process::Command; // Run programs
use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions

#[test]
fn stats() -> Result<(), Box<dyn std::error::Error>> {
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

#[test]
fn trim() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("seq_tools")?;

    cmd.arg("trim").arg("-i").arg("tests/data/reads.fastq.gz")
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