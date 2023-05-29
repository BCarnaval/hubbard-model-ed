// This is a comment, and is ignored by the compiler.
// You can test this code by clicking the "Run" button over there ->
// or if you prefer to use your keyboard, you can use the "Ctrl + Enter"
// shortcut.

use ndarray::{arr1, s, Array2};
use std::process::exit;

pub fn build_tri_up_array(diag: Vec<i32>, mut tri: Vec<i32>) -> Array2<i32> {
    // Verifying triangular number for matrix dimensions
    let dim: usize = diag.len();
    let tri_number: i32 = (dim * (dim + 1) / 2) as i32;

    if tri_number != tri.len() as i32 {
        println!(
            "Number of elements inside '{}' array does not match matrix dimension: {}.",
            String::from("tri"),
            diag.len()
        );
        exit(1)
    } else {
        // Initializing upper triangle matrix by setting the diagonal
        let mut diag_array: Array2<i32> = Array2::from_diag(&arr1(&diag));

        // Loop over upper triangle indices
        for i in 0..diag.len() {
            diag_array
                .slice_mut(s![i, (i + 1)..])
                .assign(&arr1(&tri[..dim - (i + 1)]));
            tri.drain(..dim - (i + 1));
        }
        diag_array
    }
}
