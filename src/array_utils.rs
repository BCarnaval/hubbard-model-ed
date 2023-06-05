// This module is used to define tools (functions) used on 1d arrays. For example,
// one can use the function 'build_tri_up_array(&array)' to debug the block generation
// because it provides a clean output that makes it easy to verify matrices.
use lapack::sspevd;
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
#[allow(dead_code)]
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

/// Diagonalization of upper triangular hermitian matrix using LAPACK 'sspevd'
/// Fortran implementation.
///
/// Examples
///
/// ```rust
/// let elements: Vec<f32> = vec![1., 0., 1., 0., 0., 1.];
/// let (exit_code, eig_vals): (i32, Vec<f32>) = lapack_diagonalization(elements);
/// println!("Exit code: {} and eigenvalues = {:?}", exit_code, eig_vals);
/// ```
pub fn lapack_diagonalization(lapack_ap_array: Vec<f32>) -> (i32, Vec<f32>) {
    // Matrix properties
    let mut elements: Vec<f32> = lapack_ap_array.clone();
    let array_order: i32 = get_matrix_dimension(lapack_ap_array.len()) as i32;
    let mut eigen_vals: Vec<f32> = vec![0.0; array_order as usize];
    let mut eigen_vects: Vec<f32> = Vec::with_capacity(array_order as usize);

    // Working array memory
    let lwork: i32 = 2 * array_order as i32;
    let liwork: i32 = 1;
    let mut work: Vec<f32> = Vec::with_capacity(lwork as usize);
    let mut iwork: Vec<i32> = Vec::with_capacity(liwork as usize);

    // Informative quantities
    let mut info: i32 = 0;

    unsafe {
        sspevd(
            b'N',
            b'U',
            array_order,
            &mut elements,
            &mut eigen_vals,
            &mut eigen_vects,
            1,
            &mut work,
            lwork,
            &mut iwork,
            liwork,
            &mut info,
        )
    }
    (info, eigen_vals)
}

#[cfg(test)]
mod tests {
    use ndarray::{arr2, Array2};
    use std::assert_eq;

    use crate::array_utils::{build_tri_up_array, get_matrix_dimension, lapack_diagonalization};

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

    #[test]
    fn check_lapack_sspevd() {
        let elements: Vec<f32> = vec![1., 0., 1., 0., 0., 1.];
        let output: (i32, Vec<f32>) = (0, vec![1., 1., 1.]);
        assert_eq!(output, lapack_diagonalization(elements))
    }
}
