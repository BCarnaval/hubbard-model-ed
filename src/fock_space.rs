// This is a comment, and is ignored by the compiler.
// You can test this code by clicking the "Run" button over there ->
// or if you prefer to use your keyboard, you can use the "Ctrl + Enter"
// shortcut.

use crate::file_utils::init_file_writter;
use itertools::Itertools;

#[derive(Debug)]
pub struct FockState {
    // Public attributes
    pub n_sites: i32,
    pub integer: i32,

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
    fn scalar(&self, state: i32) -> i32 {
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
        let mask: i32 = 1 << index;

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
        let mask: i32 = 1 << index;

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
    pub n_sites: i32,
    pub t: i32,
    pub u: i32,
}

impl Hubbard {
    pub fn interaction_term(&self, state_0: i32, state_1: i32) -> i32 {
        // Initializing matrix element
        let mut coefficient: i32 = 0;

        // Main loop over number of sites in the cluster (i)
        for site in 0..self.n_sites as u32 {
            // Initializing 'ket'
            let mut ket_state: FockState = FockState {
                n_sites: self.n_sites,
                integer: state_0,
                is_null: false,
            };

            // Computing 'on site' interaction using number operator
            ket_state.number(site);
            ket_state.number(site + self.n_sites as u32);

            // Updating matrix element value
            coefficient += self.u * ket_state.scalar(state_1);
        }
        coefficient
    }

    pub fn kinetic_term(&self, state_0: i32, state_1: i32) -> i32 {
        // Initializing matrix element
        let mut coefficient: i32 = 0;
        let sites: Vec<u32> = (0..self.n_sites as u32).collect();

        // Main loop over number of sites in the cluster (i, j)
        for perms in sites.iter().permutations(sites.len()).unique() {
            // Defining neighbours coordinates using permutations
            let site_i: u32 = *perms[0];
            let site_j: u32 = *perms[1];

            // Initializing 'kets' (spin up & down)
            let mut ket_state_up: FockState = FockState {
                n_sites: self.n_sites,
                integer: state_0,
                is_null: false,
            };

            let mut ket_state_down: FockState = FockState {
                n_sites: self.n_sites,
                integer: state_0,
                is_null: false,
            };

            // Kinetic term (spin up & down)
            ket_state_up.destroy(site_j);
            ket_state_up.create(site_i);

            ket_state_down.destroy(site_j + self.n_sites as u32);
            ket_state_down.create(site_i + self.n_sites as u32);

            // Updating matrix element value (spin up & down)
            coefficient += self.t * ket_state_down.scalar(state_1);
            coefficient += self.t * ket_state_up.scalar(state_1);
        }
        coefficient
    }

    pub fn get_hamiltonian(&self) {
        // Data file initialization
        let matrix_path = String::from("./Data/hubbard_hamiltonian.dat");
        let mut hubbard_wtr = init_file_writter(&matrix_path, false);

        // Main loop over Fock space states (4^(n_sites))
        for state_i in 0..(4 as u32).pow(self.n_sites as u32) {
            let mut row: Vec<i32> = Vec::new();

            for state_j in 0..(4 as u32).pow(self.n_sites as u32) {
                let mut coefficient: i32 = 0;

                // Computing 'on site' interaction
                coefficient += self.interaction_term(state_i as i32, state_j as i32);

                // Computing 'hopping/kinetic' terms
                coefficient += self.kinetic_term(state_i as i32, state_j as i32);

                // Building matrix row
                row.push(coefficient)
            }
            // Write row to data file
            hubbard_wtr.serialize(row).unwrap()
        }
    }
}
