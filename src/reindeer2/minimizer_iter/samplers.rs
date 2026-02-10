use serde::{Deserialize, Serialize};

pub trait Sampler {
    fn filter(&self, minimizer_and_kmer: (u64, u64)) -> bool;
}

#[derive(Serialize, Deserialize, Clone)]
#[cfg_attr(any(test, debug_assertions), derive(PartialEq, Debug))]
pub struct NoSampler {}

impl NoSampler {
    pub fn new() -> Self {
        Self {}
    }
}

impl Sampler for NoSampler {
    fn filter(&self, _minimizer_and_kmer: (u64, u64)) -> bool {
        true
    }
}

#[derive(Serialize, Deserialize, Clone)]
#[cfg_attr(any(test, debug_assertions), derive(PartialEq, Debug))]
pub struct MinimizerSampler {
    minimizer_mask: u64,
}

impl MinimizerSampler {
    pub fn new(last_bits_to_zero: u64) -> Self {
        let minimizer_mask = (1u64 << last_bits_to_zero) - 1;
        Self { minimizer_mask }
    }
}

impl Sampler for MinimizerSampler {
    fn filter(&self, kmer_and_minimizer: (u64, u64)) -> bool {
        dbg!(kmer_and_minimizer);
        kmer_and_minimizer.1 & self.minimizer_mask == 0
    }
}

#[derive(Serialize, Deserialize, Clone)]
#[cfg_attr(any(test, debug_assertions), derive(PartialEq, Debug))]
pub struct KmerSampler {
    kmer_mask: u64,
}

impl KmerSampler {
    pub fn new(last_bits_to_zero: u64) -> Self {
        let kmer_mask = (1u64 << last_bits_to_zero) - 1;
        Self { kmer_mask }
    }
}

impl Sampler for KmerSampler {
    fn filter(&self, kmer_and_minimizer: (u64, u64)) -> bool {
        kmer_and_minimizer.0 & self.kmer_mask == 0
    }
}
