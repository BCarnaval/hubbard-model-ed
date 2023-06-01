// This is a comment, and is ignored by the compiler.
// You can test this code by clicking the "Run" button over there ->
// or if you prefer to use your keyboard, you can use the "Ctrl + Enter"
// shortcut.

use ndarray::Array2;

/// Computes the dimension (N) of an upper triangle matrix using the 1d array length
/// of it's upper triangle elements. The equation for the length (l) is given by
/// l = N*(N + 1)/2.
///
/// Examples
///
/// ```rust
/// let ap_array: Vec<f32> = vec![1., 2., 3., 4., 5., 6., 7., 8., 9., 10.];
/// let ap_array_length: usize = ap_array.len();
/// println!("Matrix has dimension NxN, where N={}", get_matrix_dimension(ap_array_length));
/// ```
pub fn get_matrix_dimension(lapack_ap_length: usize) -> usize {
    let dimension: usize = (((1. + 8. * lapack_ap_length as f32).sqrt() - 1.) / 2.) as usize;
    dimension
}

/// Build an upper triangle matrix (A) out of the elements of it's triangle sorted using
/// LAPACK 'column-wise-packing' where A(i, j) = elements(i + j*(j - 1)/2).
///
/// Exemples
///
/// ```rust
/// let ap_array: Vec<f32> = vec![1., 2., 3., 4., 5., 6., 7., 8., 9., 10.];
/// println!("Matrix A is defined in 2d by the array:\n{:?}", build_tri_up_array(&ap_array));
/// ```
pub fn build_tri_up_array(matrix_elements: &Vec<f32>) -> Array2<f32> {
    // Verifying triangular number for matrix dimensions
    let matrix_dim: usize = get_matrix_dimension(matrix_elements.len());

    // Initializing upper triangle matrix by setting the diagonal
    let mut array: Array2<f32> = Array2::zeros((matrix_dim, matrix_dim));

    // Loop over upper triangle indices
    for i in 1..matrix_dim + 1 {
        for j in i..matrix_dim + 1 {
            let idx: usize = i + j * (j - 1) / 2;
            array[[i - 1, j - 1]] = matrix_elements[idx - 1];
        }
    }
    array
}

#[cfg(test)]
mod tests {
    use ndarray::{arr2, Array2};
    use std::assert_eq;

    use crate::array_utils::{build_tri_up_array, get_matrix_dimension};

    #[test]
    fn check_matrix_dimension() {
        let unwraped_array_length: usize = 45;
        let matrix_size: usize = 9;
        assert_eq!(matrix_size, get_matrix_dimension(unwraped_array_length));
    }

    #[test]
    fn check_tri_up_array() {
        let elements: Vec<f32> = vec![1., 2., 3., 4., 5., 6., 7., 8., 9., 10.];
        let array: Array2<f32> = arr2(&[
            [1., 2., 4., 7.],
            [0., 3., 5., 8.],
            [0., 0., 6., 9.],
            [0., 0., 0., 10.],
        ]);
        assert_eq!(array, build_tri_up_array(&elements))
    }
}
