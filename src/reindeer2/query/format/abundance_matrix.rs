use std::collections::HashMap;
use std::io::{self, Write};
use std::path::Path;

use bio::io::fasta;
use pelt::pelt;

use super::count_to_string;
use crate::reindeer2::query::{strip_header, write_header_matrix};
use crate::reindeer2::read_indexed_file_names;

use super::super::load_kmer_counts_vector;

/// Write abundances per kmer like RD1.
pub fn write_abundance_matrix(
    bf_dir: &str,
    sequence_results: &HashMap<usize, Vec<Vec<u16>>>,
    batch: &[fasta::Record],
    color_number: usize,
    breakpoints: Option<f64>,
    normalize: bool,
    writer: &mut impl Write,
) -> io::Result<()> {
    let indexed_files = Path::new(&bf_dir).join("path.txt");
    let indexed_files: Vec<String> = read_indexed_file_names(indexed_files);
    let sep = " ";
    // we need the count of kmers if we want to normalize them
    let kmer_counts = if normalize {
        load_kmer_counts_vector(bf_dir).expect("Failed to load from disk the kmer counts vector")
    } else {
        vec![color_number, 0] // TODO bizarre
    };

    write_header_matrix(writer, indexed_files, sep)?;

    for (record_id, color_vectors) in sequence_results {
        let record = &batch[*record_id];
        let id = record.id();
        let desc = record.desc().unwrap_or("");
        // Build the header string only once per record
        let seq_header = format!(">{} {}", id, desc).trim().to_string();
        let header = strip_header(&seq_header);
        write!(writer, "{header}")?;
        for (color_idx, abund_values) in color_vectors.iter().enumerate() {
            // new color => a separator
            write!(writer, "{sep}")?;

            if let Some(penalty) = breakpoints {
                let abund_values: Vec<u64> = abund_values.iter().copied().map(u64::from).collect();
                let breakpoints = pelt(&abund_values, pelt::score, penalty);
                let s = breakpoints
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",");
                write!(writer, "{}", s)?;
            } else {
                let mut start = 0;
                let mut current = abund_values[0];

                for i in 1..=abund_values.len() {
                    if i == abund_values.len() || abund_values[i] != current {
                        // new different value or end of query => we must write
                        let val_str = count_to_string(current, normalize, &kmer_counts, color_idx);
                        if start + 1 == i {
                            // only one value
                            write!(writer, "{}:{}", start, val_str)?;
                        } else {
                            // multiple values
                            write!(writer, "{}-{}:{}", start, i - 1, val_str)?;
                        }
                        if i < abund_values.len() {
                            // not the end of query
                            write!(writer, ",")?;
                            start = i;
                            current = abund_values[i];
                        }
                    }
                }
            }
        }
        writeln!(writer)?;
    }
    Ok(())
}
