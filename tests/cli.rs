
use std::process::Command; // Run programs
use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions

#[test]
fn test_stats() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("seq_tools")?;

    cmd.arg("stats").arg("-i").arg("reads.fna");
    cmd.assert()
        .stdout(predicate::str::contains("Number of nucleotides: 490"))
        .stdout(predicate::str::contains("Number of sequences: 10"))
        .stdout(predicate::str::contains("Minimum sequence length: 49"))
        .stdout(predicate::str::contains("Average sequence length: 49"));

    Ok(())
}
