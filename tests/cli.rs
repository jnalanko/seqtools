
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
