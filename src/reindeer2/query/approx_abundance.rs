fn approximate_value(log_value: u16, base: f64) -> u16 {
    if base <= 0.0 {
        panic!("base must be greater than 0");
    }
    let threshold = 1.0 / (base - 1.0);
    let logf = log_value as f64;
    if logf < threshold {
        log_value + 1
    } else {
        (base.powf((logf + 1.0) - threshold) * threshold) as u16
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
/// An approximated abundance.
///
/// Encodes the result of a query on an index, which can be:
/// - not queried (due to sampling)
/// - queried but absent
/// - queried and present
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
    const NOT_QUERIED: u16 = 0;
    const QUERIED_BUT_ABSENT: u16 = 1;

    /// Builds an approximated abundance from the result of the query on a dense index.
    ///
    /// A `dense_result` of 0 is mapped to `absent`.
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

    /// Builds an approximated abundance from the result of the query on a filter.
    ///
    /// The position of hit is the position of the 1 in the element's slice of the filter.
    /// The value is then approximated using the base.
    pub fn from_position_of_hit_in_the_filter(position_of_hit: u16, base: f64) -> Self {
        let value = approximate_value(position_of_hit, base) + 2; // FIXME limit this
        Self { value }
    }

    #[cfg(test)]
    /// Builds an approximated abundance from a raw value.
    ///
    /// This allows to acces the inner value and thus breaks the encapsulation, thus this is only allowed during tests.
    pub const fn new(val: u16) -> Self {
        Self { value: val + 2 }
    }

    /// Builds an approximated abundance from a value not queried (e.g. due to sampling).
    pub const fn new_not_queried() -> Self {
        Self {
            value: Self::NOT_QUERIED,
        }
    }

    /// Returns true if and only if the element was queried.
    pub const fn is_queried(&self) -> bool {
        self.value != Self::NOT_QUERIED
    }

    /// Returns the abundance represented.
    ///
    /// Returns
    /// - None if the element was not queried
    /// - Some(0) if the element was absent
    /// - Some(x) if the element has an approximated abundance of x
    pub const fn to_value(self) -> Option<u16> {
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

    /// Assigns an abundance to an elements given a range of possible values.
    ///
    /// Because of collisions in the filter, an element can have mutiple values.
    /// This function returns the max of these possible values, so that no collision can lead to an underestimation.
    pub fn select_abundance_from_candidates(candidates: &[(u32, Self)]) -> Option<&(u32, Self)> {
        candidates.iter().max()
    }

    /// Builds an approximated abundance from an element queried but absent.
    pub const fn new_absent() -> Self {
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
    fn test_approximate_value_with_positive_values() {
        let result = approximate_value(3, 2.0);
        assert_eq!(result, 8, "expected 2^3 to be 8");
    }

    #[test]
    #[should_panic(expected = "base must be greater than 0")]
    fn test_approximate_value_with_zero_base() {
        approximate_value(3, 0.0);
    }

    #[test]
    fn test_approximate_value_with_fractional_base() {
        let base = 2.0f64.sqrt();
        // with sqrt(2), when the approximation increase by 2, the abundance should double
        // since log_{sqrt(2)} (x) = 2 log_2(x)
        let result1 = approximate_value(6, base);
        let result2 = approximate_value(8, base);
        let result3 = approximate_value(10, base);
        assert_eq!(
            result2 / 2,
            result1,
            "check approximation consistency with base sqrt(2)"
        );
        assert_eq!(
            result3 / 2,
            result2,
            "check approximation consistency with base sqrt(2)"
        );
    }

    #[test]
    fn approx_back_and_forth() {
        for i in [0, 4, 8, 9, 4, 456, 789, 54] {
            let approx = ApproxAbundance::new(i);
            assert_eq!(approx.to_value().unwrap(), i);
        }
    }
}
