use crate::reindeer2::query::ApproxAbundance;
mod sliding_window_minimum;

pub fn fimpera(
    mut smer_abundances: Vec<Vec<Vec<ApproxAbundance>>>,
    z: usize,
) -> Vec<Vec<Vec<ApproxAbundance>>> {
    for color_vector in &mut smer_abundances {
        for smers_abondance in color_vector.iter_mut() {
            let inner = std::mem::take(smers_abondance);
            *smers_abondance = sliding_window_minimum::sliding_window_minimum(inner, z + 1)
                .expect("this should not panic");
        }
    }
    smer_abundances
}
