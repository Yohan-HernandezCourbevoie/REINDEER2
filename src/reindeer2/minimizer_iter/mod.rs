mod samplers;

use nthash::{NtHashForwardIterator, NtHashIterator};
use std::collections::VecDeque;
use thiserror::Error;

pub use samplers::{MinimizerSampler, NoSampler, Sampler};

// ///  Iterator over (k-mer, minimizer) pairs for a given sequence
// pub struct KmerMinimizerIterator<'a> {
//     seq: &'a [u8],
//     minima: Vec<u64>, // minimizers per starting position of kmers
//     current: usize,   // current k-mer start position
//     k: usize,
// }

// impl<'a> Iterator for KmerMinimizerIterator<'a> {
//     type Item = (&'a [u8], u64); // (kmer, iterator)

//     fn next(&mut self) -> Option<Self::Item> {
//         if self.current <= self.seq.len().saturating_sub(self.k) {
//             // below or equal to the last valid starting index
//             let kmer = &self.seq[self.current..self.current + self.k];
//             let minimizer = self.minima[self.current];
//             self.current += 1;
//             Some((kmer, minimizer))
//         } else {
//             None
//         }
//     }
// }

#[derive(Debug, Error)]
pub enum KmerMinimizerIteratorError {
    #[error("Sequence must be at least k (= {k}) bases long")]
    SequenceTooSmall { k: usize },
}

// returns an iterator over all the (k-mer, minimizer) from sequence input
fn kmer_minimizers_all<'a>(
    seq: &'a [u8],
    k: usize,
    m: usize,
    canonical: bool,
) -> Result<impl Iterator<Item = (u64, u64)> + 'a, KmerMinimizerIteratorError> {
    if seq.len() < k {
        return Err(KmerMinimizerIteratorError::SequenceTooSmall { k });
    }

    assert!(k >= m, "k must be greater than or equal to m");

    //  collects the hash of every m-mer in the sequence w/ rolling hash
    let m_hashes: Vec<u64> = if canonical {
        NtHashIterator::new(seq, m)
            .expect("should have been able to create canonical hash iterator")
            .collect()
    } else {
        NtHashForwardIterator::new(seq, m)
            .expect("should have been able to create hash iterator")
            .collect()
    };

    // compute sliding window to find minimums over m-mer hashes
    // for each k-mer starting at position i, the m-mers inside it are:
    //   m_hashes[i .. i + (k - m + 1)]
    let window_size = k - m + 1;
    // number of k-mers is seq.len() - k + 1; we select 1 minimizer per kmer, so there are m_hashes.len() - window_size + 1 minimizers
    let minimizers = m_hashes.len().saturating_sub(window_size) + 1;
    let mut minima = Vec::with_capacity(minimizers);
    let mut deque: VecDeque<usize> = VecDeque::new();

    // process the first window: indices from 0 .. window_size
    for i in 0..window_size {
        // positions in the window
        while let Some(&back) = deque.back() {
            // if there's a value at the back of the queue -> back is the last position P recorded in the queue
            if m_hashes[back] > m_hashes[i] {
                // if the minimizer list at position P is greater than what it is at the current position => remove from the queue
                deque.pop_back(); // remove large values from the back
            } else {
                break;
            }
        }
        deque.push_back(i); // current position i becomes a candidate for being a minimizer
    }
    // ===> the minimum for the first window has its position recorded at the front of the deque

    /* example
    window_size = 5
    m_hashes = [4, 2, 5, 1, 3]       hash values for positions 0 through 4.
    for i = 0
        deque init : [] -> deque.push_back(0) // [0]
    for i = 1
        m_hashes[0] = 4, with m_hashes[1] = 2
        4 > 2, pop index 0 from the deque
        deque.push_back(1) // [1]
    i = 2
        m_hashes[1] = 2 h m_hashes[2] = 5
        2 < 5 => break
        deque.push_back(2) // [1,2]
    i=3
        m_hashes[2] = 5 m_hashes[3] = 1
        5 > 3, pop index 2 from deque// [1]
        m_hashes[1] = 2, m_hashes[3] = 1
        2 > 1, pop index 1 // []
        deque.push_back(3) //[3]
    i=4
        m_hashes[3] = 1 m_hashes[4] = 3
        1 < 3 , break
        deque.push_back(4) //[3,4]
    => minimizer at pos 0 of deque


    */

    if let Some(&front) = deque.front() {
        minima.push(m_hashes[front]);
    }
    // then the idea carries on, we just have to remove the leftmost index that is no longer in the window each time

    //process the rest of the windows
    for i in window_size..m_hashes.len() {
        // remove indices that are now outside the current window
        while let Some(&front) = deque.front() {
            if front <= i - window_size {
                deque.pop_front();
            } else {
                break;
            }
        }
        // remove elements that are larger than the current element
        while let Some(&back) = deque.back() {
            if m_hashes[back] > m_hashes[i] {
                deque.pop_back();
            } else {
                break;
            }
        }
        deque.push_back(i);
        if let Some(&front) = deque.front() {
            minima.push(m_hashes[front]);
        }
    }

    // the number of minima should equal the number of k-mers
    debug_assert_eq!(minima.len(), seq.len() - k + 1);

    let kmer_hash_iter: Box<dyn Iterator<Item = u64>> = if canonical {
        Box::new(
            NtHashIterator::new(seq, k)
                .expect("should have been able to create canonical hash iterator"),
        )
    } else {
        Box::new(
            NtHashForwardIterator::new(seq, k)
                .expect("should have been able to create hash iterator"),
        )
    };

    // return an iterator over (hashed k-mers,corresponding minimizers)
    Ok(kmer_hash_iter.zip(minima))
}

// --- MISC ---

// fn _display_progress(total_kmers: u64, start_time: Instant) {
//     let elapsed_s = start_time.elapsed().as_secs_f64();
//     let elapsed_ms = start_time.elapsed().as_millis() as f64;
//     let kmers_per_ms = total_kmers as f64 / elapsed_ms;
//     println!(
//         "Processed: {} k-mers | Time elapsed: {:.2} s | Rate: {:.2} k-mers/ms",
//         total_kmers.to_formatted_string(&Locale::en),
//         elapsed_s,
//         kmers_per_ms
//     );
// }

/// Iterates over the (kmer, minimizer) of a sequence
pub fn kmer_minimizers_sampled<'a, S>(
    seq: &'a [u8],
    k: usize,
    m: usize,
    canonical: bool,
    sampler: &'a S,
) -> Result<impl Iterator<Item = (u64, u64)> + 'a, KmerMinimizerIteratorError>
where
    S: Sampler,
{
    let filtered = kmer_minimizers_all(seq, k, m, canonical)?
        .filter(|kmer_minimizer| sampler.filter(*kmer_minimizer));
    Ok(filtered)
}

// pub fn kmer_minimizers_sampled_marked<'a, S>(
//     seq: &'a [u8],
//     k: usize,
//     m: usize,
//     canonical: bool,
//     sampler: &'a S,
// ) -> Result<impl Iterator<Item = Option<(u64, u64)>> + 'a, KmerMinimizerIteratorError>
// where
//     S: Sampler,
// {
//     let filtered = kmer_minimizers_all(seq, k, m, canonical)?.map(|kmer_minimizer| {
//         if sampler.filter(kmer_minimizer) {
//             Some(kmer_minimizer)
//         } else {
//             None
//         }
//     });
//     Ok(filtered)
// }

#[cfg(test)]
mod tests {

    use crate::reindeer2::minimizer_iter::samplers::{KmerSampler, MinimizerSampler};

    use super::*;

    #[test]
    fn test_kmer_hash_minimizers() {
        let seq_str = "ACGTACGTACGTACGT";
        let seq_bytes = seq_str.as_bytes();

        let k = 7;
        let m = 3;
        let canonical = true;

        let no_sample_builder = NoSampler::new();

        let actual_count = kmer_minimizers_sampled(seq_bytes, k, m, canonical, &no_sample_builder)
            .unwrap()
            .count();

        // there should be seq.len() - k + 1 pairs.
        assert_eq!(actual_count, seq_bytes.len() - k + 1);
    }

    #[test]
    fn test_minimizer_sampling() {
        let seq_str = "CTTAATATCGTCCGAAAGAGTCACAGATGTAAAAGGGCGCAGTACCTAACAGGGATTGACCGACGGAACA\
        CCCTGGTCGACGCCTCAAGGCGGTAACGCCCAGTCCAATAGAGCACAATATTAGAATGCCGTCCGGCATCGCTACAGTAATTGGTGACTCG\
        TATATTACATAGGCTGTATCGACCCTGATCTTAGACGAGGTGGATTAGAAACCCAGCTCCATCGCAAGTACTAAGGTTCCGGTAACAGAAC\
        AGGCCTAGACTGAATGGTGACCATGTGCCTACCGAATTCTGGAACGTTGGGTGGCGCTTGGCTTTAGCTGCTACCATCTACCATACGGTAT\
        AACGTGGGGGTACAGCCAAAGGCAATGACCTAACGGGGCTTTTAAAGGGGGTAGATGTGCCTCGGTTGGACGGGATATGGGGGCTGTTTGT\
        CTGCCCGTGCTCATCTGCGTCTTTTTTGTATCCTTAAAATCAGTCCTGACTGCGGGGTGGCCACCGTCCACCAACGTCAGTTCTTGACTTT\
        CAAGGCCAACAGCAAGGCCTTTTGTTGGCGCTTAAGAGGAATTGGAAAGGTCTTAAACACTGCCACGGTCCACCG";
        let len = seq_str.len();
        let seq_bytes = seq_str.as_bytes();

        let k = 7;
        let m = 3;
        let canonical = true;

        let minimizer_sample_builder = KmerSampler::new(2);

        let actual_count =
            kmer_minimizers_sampled(seq_bytes, k, m, canonical, &minimizer_sample_builder)
                .unwrap()
                .count();

        // there should be seq.len() - k + 1 pairs if not sampling
        // let's check we get a quarter of that, +/- 5%
        let original_expectation = seq_bytes.len() - k + 1;
        assert!(actual_count < (original_expectation / 4) + (len * 5) / 100);
        assert!(actual_count > (original_expectation / 4) - (len * 5) / 100);
    }

    // #[test]
    // fn test_minimizer_marked() {
    //     let seq_str = "CTTAATATCGTCCGAAAGAGTCACAGATGTAAAAGGGCGCAGTACCTAACAGGGATTGACCGACGGAACA\
    //     CCCTGGTCGACGCCTCAAGGCGGTAACGCCCAGTCCAATAGAGCACAATATTAGAATGCCGTCCGGCATCGCTACAGTAATTGGTGACTCG\
    //     TATATTACATAGGCTGTATCGACCCTGATCTTAGACGAGGTGGATTAGAAACCCAGCTCCATCGCAAGTACTAAGGTTCCGGTAACAGAAC\
    //     AGGCCTAGACTGAATGGTGACCATGTGCCTACCGAATTCTGGAACGTTGGGTGGCGCTTGGCTTTAGCTGCTACCATCTACCATACGGTAT\
    //     AACGTGGGGGTACAGCCAAAGGCAATGACCTAACGGGGCTTTTAAAGGGGGTAGATGTGCCTCGGTTGGACGGGATATGGGGGCTGTTTGT\
    //     CTGCCCGTGCTCATCTGCGTCTTTTTTGTATCCTTAAAATCAGTCCTGACTGCGGGGTGGCCACCGTCCACCAACGTCAGTTCTTGACTTT\
    //     CAAGGCCAACAGCAAGGCCTTTTGTTGGCGCTTAAGAGGAATTGGAAAGGTCTTAAACACTGCCACGGTCCACCG";
    //     let seq_bytes = seq_str.as_bytes();

    //     let k = 7;
    //     let m = 3;
    //     let canonical = true;

    //     let minimizer_sample_builder = MinimizerSampler::new(2);

    //     let kmers_minis_marked =
    //         kmer_minimizers_sampled_marked(seq_bytes, k, m, canonical, &minimizer_sample_builder)
    //             .unwrap();
    //     let count = kmers_minis_marked.count();

    //     let kmers_minis_marked =
    //         kmer_minimizers_sampled_marked(seq_bytes, k, m, canonical, &minimizer_sample_builder)
    //             .unwrap();
    //     let present = kmers_minis_marked.flatten().count();

    //     let present_sampled =
    //         kmer_minimizers_sampled(seq_bytes, k, m, canonical, &minimizer_sample_builder)
    //             .unwrap()
    //             .count();

    //     // there should be seq.len() - k + 1 pairs if not sampling
    //     // let's check we get a quarter of that, +/- 5%
    //     let original_expectation = seq_bytes.len() - k + 1;

    //     assert_eq!(count, original_expectation);

    //     assert_eq!(present, present_sampled);
    // }
}
