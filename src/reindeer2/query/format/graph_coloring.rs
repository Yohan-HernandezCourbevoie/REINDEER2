use std::collections::HashMap;
use std::io::Write;

use bio::io::fasta;

use super::super::load_kmer_counts_vector;
use super::compute_median;

// rewrites a bcalm-like graph so that headers have abund info (one of the possible query operations)
pub fn graph_coloring(
    bf_dir: &str,
    sequence_results: &HashMap<usize, Vec<Vec<u16>>>,
    batch: &[fasta::Record],
    color_number: usize,
    normalize: bool,
    writer: &mut impl Write,
) -> std::io::Result<()> {
    let msg_write = "should have been able to write the query results";
    let kmer_counts = if normalize {
        load_kmer_counts_vector(bf_dir).expect("Failed to load from disk the kmer counts vector")
    } else {
        vec![color_number, 0] // TODO bizarre
    };

    for (record_id, record) in batch.iter().enumerate() {
        let id = record.id();
        let desc = record.desc().unwrap_or("");
        let full_header = format!(">{} {}", id, desc).trim().to_string();
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
            let median = if normalize {
                median as f64 / kmer_counts[color_idx] as f64 * 1_000_000f64
            } else {
                median as f64
            };
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
