// This is a comment, and is ignored by the compiler.
// You can test this code by clicking the "Run" button over there ->
// or if you prefer to use your keyboard, you can use the "Ctrl + Enter"
// shortcut.

use std::{println, vec};

use crate::array_utils::{build_tri_up_array, get_matrix_dimension};
use itertools::Itertools;
use lapack::sspevd;
use ndarray::Array2;

#[derive(Debug)]
pub struct FockState {
    // Public attributes
    pub n_sites: u32,
    pub integer: u32,

    // Private attributes
    is_null: bool,
}

impl FockState {
    /// Defines the scalar product between FockState objects.
    ///
    /// Examples
    ///
    /// ```rust
    /// let state0 = FockState { n_sites: 2, integer: 5, is_null: false };
    /// let state1 = FockState { n_sites: 2, integer: 6, is_null: false };
    /// println!("{}", state0.scalar(state1));
    /// ```
    fn scalar(&self, state: u32) -> i32 {
        // Scalar product initialization
        let mut scalar: i32 = 0;

        // Verifiying if orthogonal states
        if self.is_null {
            scalar = 0;
        } else if self.integer == state {
            scalar = 1;
        }
        scalar
    }

    /// Second quantization creation operator definition.
    ///
    /// Examples
    ///
    /// ```rust
    /// let state0 = FockState { n_sites: 2, integer: 5, is_null: false };
    /// state0.create(1);
    /// println!("{:?}", state0.integer_to_binary());
    /// ```
    fn create(&mut self, index: u32) {
        // Defining binary mask
        let mask: u32 = 1 << index;

        // Verifying if a fermion if already at position 'index' in state
        if self.integer & mask != 0 || self.is_null {
            self.is_null = true;

        // Updating Fock state integer after creating fermion
        } else {
            self.integer ^= mask;
        }
    }

    /// Second quantization anihilation operator definition.
    ///
    /// Examples
    ///
    /// ```rust
    /// let state0 = FockState { n_sites: 2, integer: 5, is_null: false };
    /// state0.destroy(1);
    /// println!("{:?}", state0.integer_to_binary());
    /// ```
    fn destroy(&mut self, index: u32) {
        // Defining binary mask
        let mask: u32 = 1 << index;

        // Verifying if no fermions are at position 'index' in state
        if self.integer & mask == 0 || self.is_null {
            self.is_null = true;

        // Updating Fock state integer after destroying fermion
        } else {
            self.integer ^= mask;
        }
    }

    /// Second quantization number operator definition.
    ///
    /// Examples
    ///
    /// ```rust
    /// let state0 = FockState { n_sites: 2, integer: 5, is_null: false };
    /// state0.number(1);
    /// println!("{:?}", state0.integer_to_binary());
    /// ```
    fn number(&mut self, index: u32) {
        // Verifiying if a fermion at site 'index' or if state is null
        if self.integer & (1 << index) == 0 || self.is_null {
            self.is_null = true;
        }
    }
}

#[derive(Debug)]
pub struct Hubbard {
    // Public attributes
    pub n_sites: u32,
    pub t: i32,
    pub u: i32,
}

impl Hubbard {
    pub fn interaction_term(&self, state_0: u32) -> i32 {
        // Initializing matrix element
        let mut coefficient: i32 = 0;

        // Main loop over number of sites in the cluster (i)
        for site in 0..self.n_sites {
            // Initializing 'ket'
            let mut ket_state: FockState = FockState {
                n_sites: self.n_sites,
                integer: state_0,
                is_null: false,
            };

            // Computing 'on site' interaction using number operator
            ket_state.number(site);
            ket_state.number(site + self.n_sites);

            // Updating matrix element value
            coefficient += self.u * ket_state.scalar(state_0);
        }
        coefficient
    }

    pub fn kinetic_term(&self, state_0: u32) -> Vec<u32> {
        // Initializing subspace states
        let mut sub_states: Vec<u32> = Vec::new();
        let sites: Vec<u32> = (0..self.n_sites).collect();

        // Main loop over number of sites in the cluster (i, j)
        for perms in sites.iter().permutations(sites.len()).unique() {
            println!("{:?}", perms);
            // Defining neighbours coordinates using permutations
            let site_i: u32 = *perms[0];
            let site_j: u32 = *perms[1];

            // Initializing 'kets' (spin up & down)
            let mut ket_up: FockState = FockState {
                n_sites: self.n_sites,
                integer: state_0,
                is_null: false,
            };

            let mut ket_down: FockState = FockState {
                n_sites: self.n_sites,
                integer: state_0,
                is_null: false,
            };

            // Kinetic term (spin up & down)
            ket_up.destroy(site_j);
            ket_up.create(site_i);

            ket_down.destroy(site_j + self.n_sites);
            ket_down.create(site_i + self.n_sites);

            // Push new state inside subspace if not already there
            if !sub_states.contains(&ket_up.integer) && !ket_up.is_null {
                sub_states.push(ket_up.integer)
            }
            if !sub_states.contains(&ket_down.integer) && !ket_down.is_null {
                sub_states.push(ket_down.integer)
            }
        }
        sub_states.sort();
        sub_states
    }

    pub fn find_sub_block(&self, state: u32) -> (Vec<u32>, Vec<f32>) {
        // Test index for new substates
        let mut idx: u32 = 0;
        let mut sub_states: Vec<u32> = vec![state];

        // Matrix elements array
        let mut elems: Vec<f32> = Vec::new();

        // Continue loop until substates aren't new
        while idx < sub_states.len() as u32 {
            // Defining current sub state
            let current_state: u32 = sub_states[idx as usize];

            // Find first hopping states
            let new_states: Vec<u32> = self.kinetic_term(current_state);
            let mut filtered: Vec<u32> = new_states.clone();
            filtered.retain(|i: &u32| !sub_states.contains(i));
            sub_states.append(&mut filtered);

            // Kinetic terms
            let states_copy: Vec<&u32> =
                sub_states.iter().filter(|i| **i < current_state).collect();
            for sub_state in states_copy {
                if !new_states.contains(sub_state) {
                    elems.push(0.);
                } else {
                    elems.push(self.t as f32);
                }
            }

            // On-site interaction coefficient
            elems.push(self.interaction_term(current_state) as f32);

            // Updating array parser
            sub_states.sort();
            idx += 1;
        }
        (sub_states, elems)
    }

    fn lapack_diagonalization(&self, lapack_ap_array: Vec<f32>) -> Vec<f32> {
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
        eigen_vals
    }

    pub fn get_eigenvalues(&self) {
        // Vector containing the blocs of the matrix & already visited states
        let mut visited: Vec<u32> = Vec::new();
        let mut blocks: Vec<Vec<u32>> = Vec::new();

        // Main loop over Fock space states (4^(n_sites))
        for state_i in 0..(4 as u32).pow(self.n_sites) {
            // Verifying if the state was already used
            if !visited.contains(&state_i) {
                // State bank from 'state_i;
                let (sub_block, matrix_elems) = self.find_sub_block(state_i);
                let eigen_vals: Vec<f32> = self.lapack_diagonalization(matrix_elems);
                println!("{:?}", eigen_vals);

                let mut filtered: Vec<u32> = sub_block.clone();

                // Building already visited states list
                filtered.retain(|i: &u32| !visited.contains(i));
                visited.append(&mut filtered);
                blocks.push(sub_block);
            } else {
                continue;
            }
        }
    }
}
