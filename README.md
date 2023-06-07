<div align="center">

# hubbard_model_ed
  
This repository contains an exact diagonalization implementation for the Hubbard model (in the approximation of 1D spins chain with first neigbohrs hopping terms using periodic boundary conditions) defined as
$$H = H_t + H_U = \sum_{i, j, \sigma}t_{ij}(c^\dagger_{i\sigma}c_{j\sigma} + h.c.) + U\sum_i n_{i\uparrow}n_{i\downarrow},$$
where $c^\dagger$ and $c$ are respectively the second quantization creation/anihilation operators and where $n$ represents the number operator from the same formalism. The code is entirely written in [Rust](https://www.rust-lang.org/) and will eventually be parallelized using [rayon](https://github.com/rayon-rs/rayon) Rust crate.

![Rust](https://img.shields.io/badge/rust-%23000000.svg?style=for-the-badge&logo=rust&logoColor=white)
  
![Build](https://img.shields.io/github/actions/workflow/status/BCarnaval/hubbard_model_ed/rust.yml?color=%23a3d1af&style=for-the-badge) ![LICENSE](https://img.shields.io/github/license/BCarnaval/hubbard_model_ed?color=blue&style=for-the-badge) ![release](https://img.shields.io/github/v/tag/BCarnaval/hubbard_model_ed?color=%23FF7F50&style=for-the-badge)
  
</div>

## Table of contents

- [Requirements](#requirements)
    - [Rust](#rust)
    - [LAPACK and BLAS](#lapack-and-blas)

- [Installation](#installation)

- [Usage](#usage)
    - [Compute eigenvalues](#compute-eigenvalues)
    - [Visualise blocks](#visualise-blocks)

- [Todo](#todo)

## Requirements

### Rust

I recommend having a version of `cargo >= 1.70.0` and [Rust](https://www.rust-lang.org/) compiler `rustc >= 1.70.0` to use this crate. If you don't know what version you have, run the following commands
```bash
$ cargo version && rustc --version
```
If you want to update them anyways, just run the command
```bash
$ rustup update
```
and it will upgrade the version of your Rust compiler, package manager, documentation and etc.

### LAPACK and BLAS

Users must also have a version of [LAPACK](https://www.netlib.org/lapack/) (Linear Algebra PACKage) and [BLAS](https://www.netlib.org/blas/) (Basic Linear Algebra Subprograms) on their computers. For Mac users, you can install it directly using [Homebrew](https://brew.sh/) and the commands
```bash
$ brew install lapack openblas
```
For additionnal details on installation I suggest to check for online support such as: [linux](https://coral.ise.lehigh.edu/jild13/2016/07/27/install-lapack-and-blas-on-linux-based-systems/), and [windows](https://icl.utk.edu/lapack-for-windows/).

## Installation

To use the program, users should clone this repository on their computer using the command
```bash
$ git clone https://github.com/BCarnaval/hubbard_model_ed ~/somewhere/hubbard_model_ed
```
Then build the binairies using `cargo` by executing the command
```bash
$ cargo build -r
```
at the root of the project. This command should use the build script to find LAPACK and BLAS and then link it inside the compiler. If the build succeed, you shoud also be able to verify if the unit tests are running properly on your machine by running the command
```bash
$ cargo test
```
at the root of the project. If all the tests pass, you are ready to use the program!

## Usage

### Compute eigenvalues

To modify the entry parameters for the 1D spins chain, users must edit the main function inside `main.rs` script. It is very simple
```rust
// Main module. This is where one should create and solve the Hubbard
// model for given parameters.
mod array_utils;
mod file_utils;
mod fock_space;

use crate::fock_space::Hubbard;
use std::println;
use std::time::Instant;

fn main() {
    let now = Instant::now();
    let hubbard_model = Hubbard {
        n_sites: 7,
        t: 1.,
        u: 2.,
    };
    hubbard_model.get_eigenvalues();
    println!("Time elapsed: {:.2?}", now.elapsed());
}
```
The parameter named `n_sites` determines how many sites are considered in the chain, the parameter `t` the hopping amplitude for the first neighbors and `u` the on-site interaction amplitude. Once the parameters are setted, run the programm using
```bash
$ cargo run -r
```
to save the eigenvalues of the hamiltonian inside `./Data/eigen_values.csv` data file.

### Visualise blocks

To visualise the different blocks of the _block diagonal hamiltonian_, one can use the function `build_tri_up_array` from the module `./src/array_utils.rs`. Simply call the function inside the module `./src/fock_space.rs` when computing matrix elements of the different blocks[^1]

[^1]: Don't forget to import the crate `array_utils.rs` at the top of `./src/fock_space.rs` in order to call the function. Use the code line: `use crate::array_utils::build_tri_up_array;`.
```rust
/// Outputs the eigenvalues of Hubbard hamiltonian by diagonalizing all
/// of it's blocks using LAPACK 'sspevd' Fortran implementation.
///
/// The eigenvalues are saved and stored inside './Data/eigen_vals.csv'.
pub fn get_eigenvalues(&self) {
    // Data file initialization (csv)
    let data_path: String = String::from("./Data/eigen_values.csv");
    let mut eig_wtr: csv::Writer<std::fs::File> = init_file_writter(&data_path, false);

    // Vector containing the blocs of the matrix & already visited states
    let mut visited: Vec<i32> = Vec::new();
    let mut blocks: Vec<Vec<i32>> = Vec::new();

    // Main loop over Fock space states (4^(n_sites))
    for state_i in 0..(4 as i32).pow(self.n_sites) {
        // Verifying if the state was already used
        if !visited.contains(&state_i) {
            // State bank from 'state_i;
            let (sub_block, matrix_elems) = self.find_sub_block(state_i);
            
            // ADD THE FOLLOWING
            println!("{:?}\n", build_tri_up_array(&matrix_elems));
            
            let (_success, eigen_vals): (i32, Vec<f32>) = lapack_diagonalization(matrix_elems);
            
            eig_wtr.serialize(eigen_vals).unwrap();

            // Building already visited states list
            let mut filtered: Vec<i32> = sub_block.clone();
            filtered.retain(|i: &i32| !visited.contains(i));
            visited.append(&mut filtered);
            blocks.push(sub_block);
        } else {
            continue;
        }
    }
}
```
and it will print all the differents blocks of the hamiltonian that are diagonalized to find the eigenvalues.

## Todo

- [x] Complete the `README.md`
- [x] Include periodic boundary conditions to hoppings operator (phase correction)
- [ ] Save eigenvalues based on filling (hamiltonian blocks)
- [ ] Save Hamiltonian blocks in text file
- [ ] Include parallel computing using [rayon](https://github.com/rayon-rs/rayon)
- [x] Comment the code base
- [x] Unit testing
