use std::cmp::{max, min};

// Local alignment of needle against the haystack.
// Returns one past the ending point of the rightmost match, if exist.
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

    for end in (1..=n).rev() {
        if score_matrix[m][end] as f64 / m as f64 >= identity_threshold {
            return Some(end);
        }
    }

    return None;
}


pub fn trim_adapters(reader: &mut impl jseqio::reader::SeqStream, output: &mut impl jseqio::writer::SeqRecordWriter, adapters: Vec<Vec<u8>>, max_trim_length: usize, min_length_after_trim: usize, identity_threshold: f64){
    while let Some(rec) = reader.read_next().unwrap() {
        let mut trim_start = 0_usize; // Trimmed read starts from there
        let mut trim_end = rec.seq.len(); // This is one past where the trimmed read ends
        for adapter in adapters.iter() {
            let start_piece = &rec.seq[0..min(max_trim_length, rec.seq.len())];
            if let Some(end) = smith_waterman(adapter, start_piece, identity_threshold) {
                trim_start = end;
            }

            let end_rev_piece: Vec<u8> = rec.seq.iter().rev().take(max_trim_length).copied().collect();
            let rev_adapter: Vec<u8> = adapter.iter().rev().copied().collect();

            if let Some(rev_end) = smith_waterman(&rev_adapter, &end_rev_piece, identity_threshold) {
                trim_end = rec.seq.len() - rev_end;
            }

        }

        let trimmed = jseqio::record::RefRecord{head: rec.head, seq: &rec.seq[trim_start..trim_end], qual: rec.qual.map(|q| &q[trim_start..trim_end])};

        if trimmed.seq.len() > min_length_after_trim {
            output.write_ref_record(&trimmed).unwrap();
        }
    }
}


#[cfg(test)]
mod tests {

    struct TestWriter {
        records: Vec<jseqio::record::OwnedRecord>
    }

    impl jseqio::writer::SeqRecordWriter for TestWriter {
        fn flush(&mut self) -> Result<(), Box<dyn std::error::Error>> {
            Ok(())  // No need to do anything
        }

        fn write_owned_record(&mut self, rec: &jseqio::record::OwnedRecord) -> Result<(), Box<dyn std::error::Error>> {
            todo!(); 
        }

        fn write_ref_record(&mut self, rec: &jseqio::record::RefRecord) -> Result<(), Box<dyn std::error::Error>> {
            self.records.push(rec.to_owned());
            Ok(())
        }
    }

    use std::io::Cursor;

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
        assert_eq!(end, 15);
    }

    #[test]
    fn test_trim_adapters(){
        let s1 =     b"TAGATACGTACGTACGTGAAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAACCGGTTAACCGGTTAACCGGTT";
        let left_adapter = b"ACTACGTACXXGT";
        let right_adapter =                                                   b"AACCGGTTAACCGGTT";
        let ans =                       b"AGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

        let mut input_fasta = Vec::<u8>::new();
        input_fasta.push(b'>');
        input_fasta.push(b'\n');
        input_fasta.extend_from_slice(s1);
        input_fasta.push(b'\n');
        
        { // Test succesful trimming
            let mut reader = jseqio::reader::DynamicFastXReader::new(Cursor::new(input_fasta.clone())).unwrap();
            let mut writer = TestWriter{records: vec![]};

            trim_adapters(&mut reader, &mut writer, vec![left_adapter.to_vec(), right_adapter.to_vec()], 50, 10, 0.65);

            assert_eq!(writer.records.len(), 1);

            eprintln!("Trimmed: {}", std::str::from_utf8(&writer.records.first().unwrap().seq).unwrap());
            assert_eq!(writer.records.first().unwrap().seq, ans);
        }

        { // Test matches too far from ends
            let mut reader = jseqio::reader::DynamicFastXReader::new(Cursor::new(input_fasta.clone())).unwrap();
            let mut writer = TestWriter{records: vec![]};
            trim_adapters(&mut reader, &mut writer, vec![left_adapter.to_vec(), right_adapter.to_vec()], 1, 10, 0.65);

            assert_eq!(writer.records.len(), 1);

            eprintln!("Trimmed: {}", std::str::from_utf8(&writer.records.first().unwrap().seq).unwrap());
            assert_eq!(writer.records.first().unwrap().seq, s1); // Unchanged from original
        }

        { // Test becomes too short
            let mut reader = jseqio::reader::DynamicFastXReader::new(Cursor::new(input_fasta)).unwrap();
            let mut writer = TestWriter{records: vec![]};
            trim_adapters(&mut reader, &mut writer, vec![left_adapter.to_vec(), right_adapter.to_vec()], 50, 100, 0.65);

            assert_eq!(writer.records.len(), 0); // Was filtered out
        }

    }
}