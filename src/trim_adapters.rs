use std::cmp::max;

// Local alignment of needle against the haystack.
// Returns one past the ending point of the leftmost match, if exist.
// Identity threshold is between 0 and 1.
fn smith_waterman(needle: &[u8], haystack: &[u8], identity_threshold: f64) -> Option<usize> {
    let m = needle.len();
    let n = haystack.len();

    let mut score_matrix = vec![vec![0_isize; n + 1]; m + 1];
    let mut max_score = 0;

    // Fill the score matrix
    for i in 1..=m {
        for j in 1..=n {
            let match_mismatch_score = if needle[i - 1] == haystack[j - 1] { 1 } else { 0 } as isize;

            score_matrix[i][j] = max(0, max(
                score_matrix[i - 1][j - 1] + match_mismatch_score,
                max(
                    score_matrix[i - 1][j] - 1,
                    score_matrix[i][j - 1] - 1,
                ),
            ));

            max_score = max(max_score, score_matrix[i][j]);
        }
    }

    // Backtrace the leftmost match
    for end in 1..=n {
        if score_matrix[m][end] as f64 / m as f64 >= identity_threshold {
            return Some(end);
        }
    }

    return None;
}


pub fn trim_adapters(reader: &mut jseqio::reader::DynamicFastXReader, output: &mut jseqio::writer::DynamicFastXWriter, adapters: Vec<Vec<u8>>, max_trim_length: usize, min_length_after_trim: usize, identity_threshold: f64){
    while let Some(rec) = reader.read_next().unwrap() {
        let ownedrec = rec.to_owned();
        for adapter in adapters.iter() {
            if let Some(end) = smith_waterman(adapter, rec.seq, identity_threshold) {
            }
        }
    }
}


#[cfg(test)]
mod tests {

    use assert_cmd::assert;

    use super::*;

    #[test]
    fn test_smith_waterman(){
        let s1 = b"TAGATACGTACGTACGTGAAG";
        let s2 =      b"ACGTAAGTACGT"; // 1 substitution
        let s3 =      b"ACTACGTACXXGT"; // See below for the optimal solution

        let end = smith_waterman(s2, s1, 0.9).unwrap();
        assert_eq!(end, 17);

        assert!(smith_waterman(s2, s1, 0.95).is_none());

        // s3 vs s1 optimal solution:
        // TAGATACGTACGTACGTGAAG;
        //      || |||||||**|*     
        //      AC-TACGTACXXGT
        // So we have 10 matches, 3 mismatches, 1 deletion -> score 10-1 = 9
        // So the largest match threshold that is positive should be 9/13.
        assert!(smith_waterman(s3, s1, 9.0/13.0 + 0.01).is_none()); // Just above the threshold
        let end = smith_waterman(s3, s1, 9.0/13.0 - 0.01).unwrap(); // Just below the threshold
        assert_eq!(end, 19);

        let s4 = b"TAC"; // Has two exact matches
        let end = smith_waterman(s4, s1, 0.999).unwrap(); // Just below the threshold
        assert_eq!(end, 7);
    }
}