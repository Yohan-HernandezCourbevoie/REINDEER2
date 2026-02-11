use crate::reindeer2::approximate_value;

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct ApproxAbundance {
    value: u16,
}

impl PartialOrd for ApproxAbundance {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ApproxAbundance {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // removing 1 by wrapping
        // so that NOT_QUERIED is mapped to `u16::MAX`
        // so when we take the min, NOT_QUERIED is not preferred
        self.value.wrapping_sub(1).cmp(&other.value.wrapping_sub(1))
    }
}

impl ApproxAbundance {
    pub const NOT_QUERIED: u16 = 0;
    pub const QUERIED_BUT_ABSENT: u16 = 1;

    pub fn from_dense(dense_result: u8, base: f64) -> Self {
        let dense_result = dense_result as u16 + 1;
        let approx_count = if dense_result == Self::NOT_QUERIED {
            Self::NOT_QUERIED
        } else if dense_result == Self::QUERIED_BUT_ABSENT {
            Self::QUERIED_BUT_ABSENT
        } else {
            approximate_value(dense_result - 2, base) + 2 // FIXME limit this
        };
        Self {
            value: approx_count,
        }
    }

    pub fn from_position_of_hit_in_the_filter(position_of_hit: u16, base: f64) -> Self {
        let value = approximate_value(position_of_hit, base) + 2; // FIXME limit this
        Self { value }
    }

    #[cfg(test)]
    pub fn new(val: u16) -> Self {
        Self { value: val + 2 }
    }

    pub fn new_not_queried() -> Self {
        Self {
            value: Self::NOT_QUERIED,
        }
    }

    pub const fn is_queried(&self) -> bool {
        self.value != Self::NOT_QUERIED
    }

    pub fn to_value(&self) -> Option<u16> {
        if self.is_queried() {
            if self.value == Self::QUERIED_BUT_ABSENT {
                Some(0)
            } else {
                Some(self.value - 2) // FIXME +1 ? -1 ?
            }
        } else {
            None
        }
    }

    pub fn select_abundance_from_candidates(candidates: &[(u32, Self)]) -> Option<&(u32, Self)> {
        candidates.iter().min()
    }

    pub fn new_absent() -> Self {
        Self {
            value: Self::QUERIED_BUT_ABSENT,
        }
    }
}

// TODO write more unit tests for this file
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn approx_back_and_forth() {
        for i in [0, 4, 8, 9, 4, 456, 789, 54] {
            let approx = ApproxAbundance::new(i);
            assert_eq!(approx.to_value().unwrap(), i);
        }
    }
}
