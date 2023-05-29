// This is a comment, and is ignored by the compiler.
// You can test this code by clicking the "Run" button over there ->
// or if you prefer to use your keyboard, you can use the "Ctrl + Enter"
// shortcut.

use std::{println, vec};

use crate::array_utils::build_tri_up_array;
use itertools::Itertools;

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
    /// std::println!("{}", state0.scalar(state1));
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
    /// std::println!("{:?}", state0.integer_to_binary());
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
    /// std::println!("{:?}", state0.integer_to_binary());
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
    /// std::println!("{:?}", state0.integer_to_binary());
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
        sub_states
    }

    fn find_sub_block(&self, state: u32) -> Vec<u32> {
        // Test index for new substates
        let mut idx: u32 = 0;
        let mut sub_states: Vec<u32> = vec![state];

        // Matrix elements array
        let mut coefficients: Vec<i32> = Vec::new();

        // Continue loop until substates aren't new
        while idx < sub_states.len() as u32 {
            // On-site interaction coefficient
            coefficients.push(self.interaction_term(sub_states[idx as usize]));

            // Find first hopping states
            let new_states: Vec<u32> = self.kinetic_term(sub_states[idx as usize]);
            let mut filtered: Vec<u32> = new_states.clone();

            // Update sub_states block by adding only new Fock states
            filtered.retain(|i: &u32| !sub_states.contains(i));
            sub_states.append(&mut filtered);

            // Add needed dimensions to resultant matrix
            if filtered.len() > 0 && idx > 0 {
                coefficients.append(&mut vec![0; filtered.len()]);
            }

            // Kinetic terms
            let cpcp = sub_states.iter().filter(|i: u32| i == state).collect();
            for sub_state in cpcp {
                if sub_state < state {
                    continue;
                } else {
                    coefficients.append(&mut vec![self.t; new_states.len()]);
                }
            }

            // Updating array parser
            idx += 1;
        }
        println!("{:?}", coefficients);
        sub_states
    }

    pub fn get_hamiltonian(&self) {
        // Vector containing the blocs of the matrix & already visited states
        let mut visited: Vec<u32> = Vec::new();
        let mut blocks: Vec<Vec<u32>> = Vec::new();

        // Main loop over Fock space states (4^(n_sites))
        for state_i in 0..(4 as u32).pow(self.n_sites) {
            // Verifying if the state was already used
            if !visited.contains(&state_i) {
                // State bank from 'state_i;
                let sub_block: Vec<u32> = self.find_sub_block(state_i);
                let mut filtered: Vec<u32> = sub_block.clone();

                // Building already visited states list
                filtered.retain(|i: &u32| !visited.contains(i));
                visited.append(&mut filtered);
                blocks.push(sub_block);
            } else {
                continue;
            }
        }
        println!("{:?}", blocks);
    }
}
