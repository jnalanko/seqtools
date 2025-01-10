use std::cmp::max;

fn backtrace_from(score_matrix: &Vec<Vec<isize>>, needle: &[u8], haystack: &[u8], mut i: usize, mut j: usize) -> usize {
    while i != 0 {
        if j == 0 { // We are at the edge. Special case to avoid accidental out of bounds access.
            assert!(score_matrix[i][j] == i as isize); // Sanity check, otherwise we have a bug
            i = 0; // The backtrace will just decrement i until done.
            break;
        }

        // Here both i > 0 and j > 0, so we all these accesses are in bounds.
        if needle[i-1] == haystack[j-1] && score_matrix[i-1][j-1] == score_matrix[i][j] - 1 { // Prefer matches
            i -= 1;
            j -= 1;
        } else if score_matrix[i-1][j-1] == score_matrix[i][j] { // Prefer diagonal moves
            i -= 1;
            j -= 1;
        } else if score_matrix[i-1][j] == score_matrix[i][j] - 1 { // Prefer vertical moves
            i -= 1;
        } else if score_matrix[i][j-1] == score_matrix[i][j] - 1 {
            j -= 1;
        } else {
            panic!("Bug in Smith-Waterman backtrace");
        }
    }
    j
}

// Local alignment of needle against the haystack.
// Returns the starting point of the leftmost and rightmost math, if exist.
// Identity threshold is between 0 and 1.
fn smith_waterman(needle: &[u8], haystack: &[u8], identity_threshold: f64) -> Option<(usize, usize)> {
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

    for line in score_matrix.iter() { // DEBUG
        eprintln!("{:?}", line);
    }

    let mut leftmost: Option<usize> = None;
    let mut rightmost: Option<usize> = None;

    // Backtrace the leftmost match
    for end in 1..=n {
        if score_matrix[m][end] as f64 / m as f64 >= identity_threshold {
            leftmost = Some(backtrace_from(&score_matrix, needle, haystack, m, end));
            break;
        }
    }

    // Backtrace the rightmost match
    for end in (1..=n).rev() {
        if score_matrix[m][end] as f64 / m as f64 >= identity_threshold {
            rightmost = Some(backtrace_from(&score_matrix, needle, haystack, m, end));
            break;
        }
    }

    if leftmost.is_some() && rightmost.is_some() {
        Some((leftmost.unwrap(), rightmost.unwrap()))
    } else {
        None
    }

}


pub fn trim_adapters(reader: &mut jseqio::reader::DynamicFastXReader, output: &mut jseqio::writer::DynamicFastXWriter, adapters: Vec<Vec<u8>>, max_length_to_seq: usize, min_length_after_trim: usize){
    todo!();
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


        let (left,right) = smith_waterman(s2, s1, 0.9).unwrap();
        assert_eq!(left, 5);
        assert_eq!(right, 5);

        assert!(smith_waterman(s2, s1, 0.95).is_none());

        // s3 vs s1 optimal solution:
        // TAGATACGTACGTACGTGAAG;
        //      || |||||||**|*     
        //      AC-TACGTACXXGT
        // So we have 10 matches, 3 mismatches, 1 deletion -> score 10-1 = 9
        // So the largest match threshold that is positive should be 9/13.
        assert!(smith_waterman(s3, s1, 9.0/13.0 + 0.01).is_none()); // Just above the threshold
        let (left, right) = smith_waterman(s3, s1, 9.0/13.0 - 0.01).unwrap(); // Just below the threshold
        assert_eq!(left, 5);
        assert_eq!(right, 5);

        // Backtrace needs to check that there really exists a match!
    }
}