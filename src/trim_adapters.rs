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

pub struct TrimStats {
    pub bases_trimmed_from_start: usize,
    pub bases_trimmed_from_end: usize,
    pub reads_with_adapter_at_start: usize,
    pub reads_with_adapter_at_end: usize,
    pub reads_with_adapter_at_both_ends: usize,
    pub discarded_reads: usize,
    pub bases_in_discarded_reads: usize,
    pub start_found_counts: Vec<usize>, // For each adapter
    pub end_found_counts: Vec<usize>,
    pub total_start_distance: Vec<usize>,
    pub total_end_distance: Vec<usize>,
}

pub fn trim_adapters(reader: &mut impl jseqio::reader::SeqStream, output: &mut impl jseqio::writer::SeqRecordWriter, adapters: Vec<Vec<u8>>, max_trim_length: usize, min_length_after_trim: usize, identity_threshold: f64) -> TrimStats {

    let mut stats = TrimStats{
        bases_trimmed_from_start: 0, 
        bases_trimmed_from_end: 0, 
        reads_with_adapter_at_start: 0, 
        reads_with_adapter_at_end: 0, 
        reads_with_adapter_at_both_ends: 0, 
        discarded_reads: 0, 
        bases_in_discarded_reads: 0, 
        start_found_counts: vec![0; adapters.len()], 
        end_found_counts: vec![0; adapters.len()], 
        total_start_distance: vec![0; adapters.len()], 
        total_end_distance: vec![0; adapters.len()]};

    let mut n_reads = 0_usize;
    let mut total_input_length = 0_usize;
    let mut total_output_length = 0_usize;

    while let Some(rec) = reader.read_next().unwrap() {
        n_reads += 1;
        total_input_length += rec.seq.len();
        let mut trim_start = 0_usize; // Trimmed read starts from there
        let mut trim_end = rec.seq.len(); // This is one past where the trimmed read ends
        for (adapter_idx, adapter) in adapters.iter().enumerate() {
            let start_piece = &rec.seq[0..min(max_trim_length, rec.seq.len())];
            if let Some(end) = smith_waterman(adapter, start_piece, identity_threshold) {
                stats.start_found_counts[adapter_idx] += 1;
                stats.total_start_distance[adapter_idx] += end;
                trim_start = end;
            }

            let end_rev_piece: Vec<u8> = rec.seq.iter().rev().take(max_trim_length).copied().collect();
            let rev_adapter: Vec<u8> = adapter.iter().rev().copied().collect();

            if let Some(rev_end) = smith_waterman(&rev_adapter, &end_rev_piece, identity_threshold) {
                stats.end_found_counts[adapter_idx] += 1;
                stats.total_end_distance[adapter_idx] += rev_end;
                trim_end = rec.seq.len() - rev_end;
            }

        }
        
        if trim_start != 0 {
            stats.reads_with_adapter_at_start += 1;
        }
        if trim_end != rec.seq.len() {
            stats.reads_with_adapter_at_end += 1;
        }
        if trim_start != 0 && trim_end != rec.seq.len() {
            stats.reads_with_adapter_at_both_ends += 1;
        }

        if trim_start > trim_end { // Overlapping trims
            trim_start = trim_end; // Remove the overlap. Everything will be trimmed.
        }

        let trimmed = jseqio::record::RefRecord{head: rec.head, seq: &rec.seq[trim_start..trim_end], qual: rec.qual.map(|q| &q[trim_start..trim_end])};

        if trimmed.seq.len() > min_length_after_trim {
            output.write_ref_record(&trimmed).unwrap();
            total_output_length += trimmed.seq.len();
            stats.bases_trimmed_from_start += trim_start;
            stats.bases_trimmed_from_end += rec.seq.len() - trim_end;
        } else {
            stats.discarded_reads += 1;
            stats.bases_in_discarded_reads += rec.seq.len();
        }
    }

    println!("Adapter\tFound near start\tFound near end\tMean distance to start\tMean distance from end");
    for (adapter_idx, adapter) in adapters.iter().enumerate() {
        println!("{}\t{}\t{}\t{}\t{}", std::str::from_utf8(adapter).unwrap(), stats.start_found_counts[adapter_idx], stats.end_found_counts[adapter_idx], stats.total_start_distance[adapter_idx] as f64 / stats.start_found_counts[adapter_idx] as f64, stats.total_end_distance[adapter_idx] as f64 / stats.end_found_counts[adapter_idx] as f64);
    }

    println!("Total number of bases in input: {}", total_input_length);
    println!("Total number of bases in output: {}", total_output_length);
    println!("Total fraction of bases removed: {:.2}%", (total_input_length - total_output_length) as f64 / total_input_length as f64 * 100.0);

    println!("Reads with adapter near start: {} ({:.2}%)", stats.reads_with_adapter_at_start, stats.reads_with_adapter_at_start as f64 / n_reads as f64 * 100.0);
    println!("Reads with adapter near end: {} ({:.2}%)", stats.reads_with_adapter_at_end, stats.reads_with_adapter_at_end as f64 / n_reads as f64 * 100.0);
    println!("Reads with adapter near both ends: {} ({:.2}%)", stats.reads_with_adapter_at_both_ends, stats.reads_with_adapter_at_both_ends as f64 / n_reads as f64 * 100.0);

    println!("Bases trimmed from starts: {}", stats.bases_trimmed_from_start);
    println!("Bases trimmed from ends: {}", stats.bases_trimmed_from_end);
    println!("Discarded reads (too short after possible trimmming): {} ({:.2}%)", stats.discarded_reads, stats.discarded_reads as f64 / n_reads as f64 * 100.0);
    println!("Bases in discarded reads: {}", stats.bases_in_discarded_reads);

    stats

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