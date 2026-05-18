mod samplers;
mod sequence_utils;

use itertools::Itertools;
pub use samplers::{KmerSampler, MinimizerSampler, NoSampler, Sampler};
use sequence_utils::{SKInfos, is_canonical, new_sk_iterator, reverse_complement};

use thiserror::Error;
use xxhash_rust::xxh3;

pub fn get_hash<const CANONICAL: bool>(seq: &[u8]) -> u64 {
    if CANONICAL {
        if is_canonical(seq) {
            xxh3::xxh3_64(seq)
        } else {
            let revcomp = reverse_complement(seq).collect_vec();
            xxh3::xxh3_64(&revcomp)
        }
    } else {
        xxh3::xxh3_64(seq)
    }
}

pub fn get_all_minimizers_and_kmers<const CANONICAL: bool, T: Iterator<Item = SKInfos>>(
    superkmers: T,
    seq: &[u8],
    k: usize,
    m: usize,
) -> impl Iterator<Item = (u64, u64)> + use<'_, CANONICAL, T> {
    superkmers.flat_map(move |sk| {
        let minimizer = &seq[sk.minimizer_start..sk.minimizer_start + m];
        let minimizer = get_hash::<CANONICAL>(minimizer);

        (sk.superkmer_start..(sk.superkmer_end - k + 1)).map(move |kmer_start| {
            let kmer = &seq[kmer_start..kmer_start + k];
            let kmer = get_hash::<CANONICAL>(kmer);
            (kmer, minimizer)
        })
    })
}

#[derive(Debug, Error)]
pub enum KmerMinimizerIteratorError {
    #[error("Sequence must be at least k (= {k}) bases long")]
    SequenceTooSmall { k: usize },
}

// returns an iterator over all the (k-mer, minimizer) from sequence input
pub fn kmer_minimizers_all<'a>(
    seq: &'a [u8],
    k: usize,
    m: usize,
    canonical: bool,
) -> Result<impl Iterator<Item = (u64, u64)> + 'a, KmerMinimizerIteratorError> {
    if seq.len() < k {
        return Err(KmerMinimizerIteratorError::SequenceTooSmall { k });
    }
    let width = (k - m + 1) as u16;
    let kmer_hash_iter: Box<dyn Iterator<Item = (u64, u64)>> = if canonical {
        let iterator = new_sk_iterator::<true>(seq, m, width);
        Box::new(get_all_minimizers_and_kmers::<true, _>(iterator, seq, k, m))
    } else {
        let iterator = new_sk_iterator::<false>(seq, m, width);
        Box::new(get_all_minimizers_and_kmers::<false, _>(
            iterator, seq, k, m,
        ))
    };
    Ok(kmer_hash_iter)
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

// OPTIMIZE we might be more efficient if we filtered the minimizer at the superkmer level
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

    use rand::seq::IndexedRandom;

    pub fn random_dna_read(len: usize) -> Vec<u8> {
        const ALPHABET: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let mut rng = rand::rng();

        let mut read = Vec::with_capacity(len);
        for _ in 0..len {
            read.push(*ALPHABET.choose(&mut rng).unwrap());
        }
        read
    }

    #[test]
    fn test_kmer_sampling() {
        let seq = random_dna_read(3000000);

        let k = 31;
        let m = 15;
        let canonical = true;
        for i in 0..10 {
            let kmer_sample_builder = KmerSampler::new(i);

            let actual_count = kmer_minimizers_sampled(&seq, k, m, canonical, &kmer_sample_builder)
                .unwrap()
                .count();

            // there should be seq.len() - k + 1 pairs if not sampling
            // let's check we get a quarter of that, +/- 5%
            let original_expectation = seq.len() - k + 1;
            let expectation = original_expectation / 2usize.pow(i as u32);

            assert!(actual_count < expectation + (expectation * 5) / 100);
            assert!(actual_count > expectation - (expectation * 5) / 100);
        }
    }

    #[test]
    fn test_minimizer_sampling() {
        let seq = random_dna_read(3000000);

        let k = 31;
        let m = 15;
        let canonical = true;
        for i in 1..10 {
            let minimizer_sampler = MinimizerSampler::new(i);

            let actual_count = kmer_minimizers_sampled(&seq, k, m, canonical, &minimizer_sampler)
                .unwrap()
                .count();

            // there should be seq.len() - k + 1 pairs if not sampling
            // let's check we get a fraction of that, +/- 5%
            let original_expectation = seq.len() - k + 1;
            let expectation = original_expectation / 2usize.pow(i as u32);

            assert!(actual_count < expectation + (expectation * 20) / 100);
            assert!(actual_count > expectation - (expectation * 20) / 100);
        }
    }

    #[test]
    fn test_get_hash() {
        let seq = "TCGAGCTCCTAACGCGCAGGTTTTTTCCCGTAACGTACGTGAGAGTATCATACACACGCAAACGGAAGTGAGCGTCAGTGTTCGAACTCGCTCTCATCTG";
        let revcomp = "CAGATGAGAGCGAGTTCGAACACTGACGCTCACTTCCGTTTGCGTGTGTATGATACTCTCACGTACGTTACGGGAAAAAACCTGCGCGTTAGGAGCTCGA";
        let seq = seq.as_bytes();
        let revcomp = revcomp.as_bytes();

        assert_eq!(get_hash::<true>(seq), get_hash::<true>(revcomp));
    }
}
