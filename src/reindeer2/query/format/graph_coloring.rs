use std::collections::HashMap;
use std::io::Write;

use bio::io::fasta;

use super::{
    compute_median, count_to_string_witout_star, get_full_header, KmerCountsAndNormalizeValue,
};

// rewrites a bcalm-like graph so that headers have abund info (one of the possible query operations)
pub fn graph_coloring(
    sequence_results: &HashMap<usize, Vec<Vec<u16>>>,
    batch: &[fasta::Record],
    normalize: &Option<KmerCountsAndNormalizeValue>,
    writer: &mut impl Write,
) -> std::io::Result<()> {
    let msg_write = "should have been able to write the query results";

    for (record_id, record) in batch.iter().enumerate() {
        let full_header = get_full_header(record);
        let seq_str = std::str::from_utf8(record.seq()).expect("Invalid UTF-8 sequence");

        // Let's load the results
        // color_vectors is a Vec<Vec<u16>>. Each index = a color,
        // each inner Vec<u16> = all abundance values for that color
        let color_vectors = sequence_results
            .get(&record_id)
            .expect("should have been able to get the result from the record id");

        // if no data, just write the original header
        if color_vectors.iter().all(Vec::is_empty) {
            writeln!(writer, "{}\n{}", full_header, seq_str).expect(msg_write);
            continue;
        }
        // otherwise, build an augmented header
        let mut header_parts = Vec::with_capacity(color_vectors.len() + 1);
        header_parts.push(full_header);

        // for each color, we do the median of all values:
        for (color_idx, vals) in color_vectors.iter().enumerate() {
            if vals.is_empty() {
                // skip color if it has no data
                continue;
            }
            let median = compute_median(vals);
            let median = count_to_string_witout_star(median, normalize, color_idx);
            // push e.g. "col:1:12"
            header_parts.push(format!("col:{}:{}", color_idx, median));
        }

        // join info like
        // ">seq1 col:0:12 col:1:29"
        let new_header = header_parts.join(" ");
        writeln!(writer, "{}\n{}", new_header, seq_str).expect(msg_write);
    }
    Ok(())
}
