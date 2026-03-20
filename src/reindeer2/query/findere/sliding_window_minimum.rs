pub fn sliding_window_minimum<T: Copy + Ord>(mut arr: Vec<T>, w: usize) -> Result<Vec<T>, Vec<T>> {
    // basically, the minimum in the sliding window is the minimum of two "offseted" vectors (min_left and min_right)
    // but:
    // min_left is computed on the fly (we only store its required element)
    // min_right is stored in the vector passed as parameter
    // the response is computed in the arr vector directly to avoid memory allocation.
    // OPTIMIZE throw errors instead of returning empty response
    let n = arr.len();

    if n < w {
        arr.truncate(0);
        return Err(arr);
    }

    if w <= 1 {
        return Ok(arr);
    }

    let nb_win = n / w;
    let nb_elem_last_window = n % w;

    let mut min_left = arr[0];
    for x in arr.iter().take(w) {
        min_left = min_left.min(*x);
    }

    for i in 0..nb_win - 1 {
        let start_window = i * w;

        // compute min_right in-place
        for indice in (start_window..=(start_window + w - 2)).rev() {
            // we compute "min_right" here, directly in arr vector
            arr[indice] = arr[indice].min(arr[indice + 1]);
        }

        for j in 0..w {
            arr[start_window + j] = arr[start_window + j].min(min_left);
            min_left = if j == 0 {
                arr[start_window + w + j]
            } else {
                min_left.min(arr[start_window + w + j])
            };
        }
    }

    // last window
    let start_window = (nb_win - 1) * w;

    // compute min_right for last window
    for indice in (start_window..=start_window + w - 2).rev() {
        arr[indice] = arr[indice].min(arr[indice + 1]);
    }

    // compute the min for the last window
    for j in 0..nb_elem_last_window {
        arr[start_window + j] = arr[start_window + j].min(min_left);
        min_left = if j == 0 {
            arr[start_window + w + j]
        } else {
            min_left.min(arr[start_window + w + j])
        };
    }

    arr[start_window + nb_elem_last_window] = arr[start_window + nb_elem_last_window].min(min_left);

    arr.truncate(arr.len() - (w - 1));
    Ok(arr)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper to run the function cleanly without mutating the original
    fn swm(arr: &[i32], w: usize) -> Vec<i32> {
        let v = arr.to_vec();
        match sliding_window_minimum(v, w) {
            Ok(v) => v,
            Err(v) => v,
        }
    }

    // --- Edge cases ---

    #[test]
    fn test_empty_array() {
        assert_eq!(swm(&[], 3), Vec::<i32>::new());
    }

    #[test]
    fn test_w_larger_than_array() {
        assert_eq!(swm(&[1, 2, 3], 5), Vec::<i32>::new());
    }

    #[test]
    fn test_w_equals_array_size() {
        assert_eq!(swm(&[3, 1, 2], 3), vec![1]);
    }

    #[test]
    fn test_w_zero() {
        assert_eq!(swm(&[1, 2, 3], 0), vec![1, 2, 3]);
    }

    #[test]
    fn test_w_one() {
        assert_eq!(swm(&[4, 2, 7, 1], 1), vec![4, 2, 7, 1]);
    }

    // --- Correctness against naive implementation ---

    fn naive_sliding_window_minimum(arr: &[i32], w: usize) -> Vec<i32> {
        if arr.len() < w || w == 0 {
            return vec![];
        }
        (0..=arr.len() - w)
            .map(|i| *arr[i..i + w].iter().min().unwrap())
            .collect()
    }

    fn check(arr: &[i32], w: usize) {
        assert_eq!(
            swm(arr, w),
            naive_sliding_window_minimum(arr, w),
            "Failed for arr={:?}, w={}",
            arr,
            w
        );
    }

    #[test]
    fn test_increasing() {
        check(&[1, 2, 3, 4, 5, 6], 3);
    }

    #[test]
    fn test_decreasing() {
        check(&[6, 5, 4, 3, 2, 1], 3);
    }

    #[test]
    fn test_all_equal() {
        check(&[3, 3, 3, 3, 3], 3);
    }

    #[test]
    fn test_single_minimum_at_start() {
        check(&[1, 5, 5, 5, 5, 5], 3);
    }

    #[test]
    fn test_single_minimum_at_end() {
        check(&[5, 5, 5, 5, 5, 1], 3);
    }

    #[test]
    fn test_single_minimum_in_middle() {
        check(&[5, 5, 1, 5, 5, 5], 3);
    }

    #[test]
    fn test_negative_values() {
        check(&[-3, -1, -4, -1, -5, -9, -2, -6], 3);
    }

    #[test]
    fn test_mixed_values() {
        check(&[3, -1, 2, 4, -2, 7, 1], 3);
    }

    // --- Block boundary stress ---

    #[test]
    fn test_exact_two_blocks() {
        check(&[3, 1, 4, 1, 5, 9], 3); // n = 2*w exactly
    }

    #[test]
    fn test_exact_three_blocks() {
        check(&[3, 1, 4, 1, 5, 9, 2, 6, 5], 3); // n = 3*w exactly
    }

    #[test]
    fn test_last_block_partial() {
        check(&[3, 1, 4, 1, 5, 9, 2], 3); // n % w == 1
    }

    #[test]
    fn test_last_block_partial_2() {
        check(&[3, 1, 4, 1, 5, 9, 2, 6], 3); // n % w == 2
    }

    #[test]
    fn test_w_two() {
        check(&[5, 3, 1, 4, 2, 6], 2);
    }

    #[test]
    fn test_large_w() {
        check(&[5, 3, 1, 4, 2, 6, 7, 8, 9, 0], 7);
    }

    // --- Randomized fuzzing against naive ---
    #[test]
    fn test_fuzz() {
        // Deterministic "random" sequence to keep tests reproducible
        let arr: Vec<i32> = (0..50).map(|i| ((i * 37 + 13) % 20) - 10).collect();
        for w in 1..=10 {
            check(&arr, w);
        }
    }
}
