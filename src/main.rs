/// Main module. This is where one should create and solve the Hubbard
/// model for given parameters.
mod array_utils;
mod file_utils;
mod fock_space;

use crate::fock_space::Hubbard;
use std::println;
use std::time::Instant;

fn main() {
    let now = Instant::now();
    let hubbard_model = Hubbard {
        n_sites: 3,
        t: 1.,
        u: 2.,
    };
    hubbard_model.find_sub_block(6);
    println!("Time elapsed: {:.2?}", now.elapsed());
}
