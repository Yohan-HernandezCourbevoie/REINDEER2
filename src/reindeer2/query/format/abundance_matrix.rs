use std::collections::HashMap;
use std::io::{self, Write};

use bio::io::fasta;
use pelt::pelt;

use super::count_to_string_with_star;
use super::KmerCountsAndNormalizeValue;
use crate::reindeer2::query::format::get_full_header;
use crate::reindeer2::query::strip_header;

pub const SEPARATOR: &str = " ";

/// Writes the header of the matrix, including a `\n` at the end
pub fn write_header_matrix(
    writer: &mut impl Write,
    indexed_files: Vec<String>,
    sep: &str,
) -> io::Result<()> {
    write!(writer, "query")?;
    for indexed_file in indexed_files {
        write!(writer, "{sep}{indexed_file}")?;
    }
    writeln!(writer)?;
    Ok(())
}

/// Write abundances per kmer like RD1.
pub fn write_abundance_matrix(
    sequence_results: &HashMap<usize, Vec<Vec<u16>>>,
    batch: &[fasta::Record],
    breakpoints: &Option<f64>,
    normalize: &Option<KmerCountsAndNormalizeValue>,
    writer: &mut impl Write,
) -> io::Result<()> {
    // we need the count of kmers if we want to normalize them

    for (record_id, color_vectors) in sequence_results {
        let record = &batch[*record_id];
        // Build the header string only once per record
        let seq_header = get_full_header(record);
        // OPTIMIZE unecessary computation here
        let header = strip_header(&seq_header);
        write!(writer, "{header}")?;
        for (color_idx, abund_values) in color_vectors.iter().enumerate() {
            // new color => a separator
            write!(writer, "{SEPARATOR}")?;

            if let Some(penalty) = breakpoints {
                let abund_values: Vec<u64> = abund_values.iter().copied().map(u64::from).collect();
                let breakpoints = pelt(&abund_values, pelt::score, *penalty);
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
                        let val_str = count_to_string_with_star(current, normalize, color_idx);
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
