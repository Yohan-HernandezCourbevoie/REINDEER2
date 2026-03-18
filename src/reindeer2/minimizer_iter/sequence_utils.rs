use simd_minimizers::canonical_minimizers;
use simd_minimizers::minimizers;
use simd_minimizers::packed_seq::PackedSeqVec;
use simd_minimizers::packed_seq::SeqVec;

use std::cmp::Ordering;
use std::fmt::Debug;
use std::iter::Copied;
use std::iter::Map;
use std::iter::Rev;
use std::slice::Iter;

// Branch prediction hint. This is currently only available on nightly but it
// consistently improves performance by 10-15%.
#[cfg(not(feature = "nightly"))]
use core::convert::identity as unlikely;
#[cfg(feature = "nightly")]
use core::intrinsics::unlikely;

pub const REVCOMP_TAB: [u8; 255] = {
    let mut tab = [0; 255];
    tab[b'A' as usize] = b'T';
    tab[b'T' as usize] = b'A';
    tab[b'C' as usize] = b'G';
    tab[b'G' as usize] = b'C';
    tab[b'N' as usize] = b'T'; // N is A, so its complement is T
    tab
};

#[allow(clippy::type_complexity)]
pub fn reverse_complement<'a>(seq: &'a [u8]) -> Map<Rev<Iter<'a, u8>>, fn(&'a u8) -> u8> {
    seq.iter()
        .rev()
        .map(|base| unsafe { *REVCOMP_TAB.get_unchecked(*base as usize) })
}

pub fn same_orientation(seq: &[u8]) -> Copied<Iter<'_, u8>> {
    seq.iter().copied()
}

pub fn is_canonical(seq: &[u8]) -> bool {
    let mut orientation_1 = same_orientation(seq);
    let mut orientation_2 = reverse_complement(seq);
    while let (Some(xc), Some(yc)) = (orientation_1.next(), orientation_2.next()) {
        let xc = if unlikely(xc == b'N') { b'A' } else { xc };

        match xc.cmp(&yc) {
            Ordering::Less => return true,
            Ordering::Greater => return false,
            Ordering::Equal => {}
        }
    }
    // in case of palindrome, prefer saying the sequence is canonical
    true
}

#[derive(Copy, Clone, Debug)]
pub struct SKInfos {
    pub minimizer_start: usize,
    pub superkmer_start: usize,
    pub superkmer_end: usize,
}

pub struct SKInfosIterator {
    pos: usize,
    k: usize,
    size_read: usize,
    min_pos_vec: Vec<u32>,
    sks_pos_vec: Vec<u32>,
}

impl SKInfosIterator {
    fn new(
        k: usize,
        size_read: usize,
        min_pos_vec: Vec<u32>,
        sks_pos_vec: Vec<u32>,
    ) -> Option<Self> {
        if min_pos_vec.len() != sks_pos_vec.len() {
            return None;
        }

        Some(Self {
            k,
            size_read,
            pos: 0,
            min_pos_vec,
            sks_pos_vec,
        })
    }
}

impl Iterator for SKInfosIterator {
    type Item = SKInfos;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.min_pos_vec.len() {
            return None;
        }

        let infos = if self.pos + 1 < self.min_pos_vec.len() {
            SKInfos {
                minimizer_start: self.min_pos_vec[self.pos] as usize,
                superkmer_start: self.sks_pos_vec[self.pos] as usize,
                superkmer_end: self.sks_pos_vec[self.pos + 1] as usize + self.k - 1,
            }
        } else {
            SKInfos {
                minimizer_start: self.min_pos_vec[self.pos] as usize,
                superkmer_start: self.sks_pos_vec[self.pos] as usize,
                superkmer_end: self.size_read,
            }
        };

        self.pos += 1;
        Some(infos)
    }
}

fn canonical_minimizer_and_superkmer_positions(
    seq: &[u8],
    k: usize,
    m: usize,
) -> (Vec<u32>, Vec<u32>) {
    let mut minimizer_positions = Vec::new();
    let mut super_kmers = Vec::new();
    let packed_seq = PackedSeqVec::from_ascii(seq);
    canonical_minimizers(m, k - m + 1)
        .super_kmers(&mut super_kmers)
        .run(packed_seq.as_slice(), &mut minimizer_positions);
    (minimizer_positions, super_kmers)
}

fn minimizer_and_superkmer_positions(seq: &[u8], k: usize, m: usize) -> (Vec<u32>, Vec<u32>) {
    let mut minimizer_positions = Vec::new();
    let mut super_kmers = Vec::new();
    let packed_seq = PackedSeqVec::from_ascii(seq);
    minimizers(m, k - m + 1)
        .super_kmers(&mut super_kmers)
        .run(packed_seq.as_slice(), &mut minimizer_positions);
    (minimizer_positions, super_kmers)
}

pub fn new_sk_iterator<const CANONICAL: bool>(
    seq: &[u8],
    minimizer_size: usize,
    width: u16,
) -> SKInfosIterator {
    assert_eq!(
        width % 2,
        1,
        "width must be odd to break ties between multiple minimizers"
    );
    let k = minimizer_size + width as usize - 1;

    let (min_pos_vec, sks_pos_vec) = if CANONICAL {
        canonical_minimizer_and_superkmer_positions(
            seq,
            k, // TODO ask
            minimizer_size,
        )
    } else {
        minimizer_and_superkmer_positions(seq, k, minimizer_size)
    };

    SKInfosIterator::new(k, seq.len(), min_pos_vec, sks_pos_vec)
        .expect("cannot build a super-k-mer over the input sequence (is the sequence too short ?)")
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    #[test]
    fn test_simd_iterator() {
        let seq = "TATACTCTTAGGGTTGGTTGCATGGTGAAGT";
        let k = seq.len();
        let m = 15;

        let minimizer = {
            let (min_pos_vec, _) =
                canonical_minimizer_and_superkmer_positions(seq.as_bytes(), k, m);
            assert_eq!(min_pos_vec.len(), 1);
            &seq[min_pos_vec[0] as usize..min_pos_vec[0] as usize + m]
        };
        let rc = "ACTTCACCATGCAACCAACCCTAAGAGTATA";

        let minimizer_rc = {
            let (min_pos_vec, _) = canonical_minimizer_and_superkmer_positions(rc.as_bytes(), k, m);
            assert_eq!(min_pos_vec.len(), 1);
            &rc[min_pos_vec[0] as usize..min_pos_vec[0] as usize + m]
        };

        assert_eq!(
            minimizer.as_bytes(),
            &reverse_complement(minimizer_rc.as_bytes()).collect_vec()
        )
    }
}
