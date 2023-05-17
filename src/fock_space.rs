// This is a comment, and is ignored by the compiler.
// You can test this code by clicking the "Run" button over there ->
// or if you prefer to use your keyboard, you can use the "Ctrl + Enter"
// shortcut.

use crate::file_utils::init_file_writter;

#[derive(Debug)]
pub struct FockState {
    // Public attributes
    pub n_sites: i32,
    pub integer: i32,

    // Private attributes
    is_null: bool,
}

impl FockState {
    /// Converts the Fock state integer to vector containing
    /// associated bits.
    ///
    /// Examples
    ///
    /// ```rust
    /// let state0 = FockState { n_sites: 2, integer: 5, is_null: false };
    /// std::println!("{:?}", state0.integer_to_binary);
    /// ```
    pub fn integer_to_binary(&self) -> Vec<i32> {
        // Initializing divided part and bits array
        let mut divided: i32 = self.integer;
        let configs: u32 = (2 as u32).pow(self.n_sites as u32);
        let mut binary_vec: Vec<i32> = vec![0; configs as usize];

        // Lopp while the 'divided' is different than 0
        let mut idx: usize = 0;
        while divided > 0 {
            let remainder: i32 = divided % 2;
            binary_vec[idx] = remainder;

            // Updating the divided part using previous value
            divided /= 2;
            idx += 1;
        }
        binary_vec
    }

    /// Defines the scalar product between FockState objects.
    ///
    /// Examples
    ///
    /// ```rust
    /// let state0 = FockState { n_sites: 2, integer: 5, is_null: false };
    /// let state1 = FockState { n_sites: 2, integer: 6, is_null: false };
    /// std::println!("{}", state0.scalar(state1));
    /// ```
    fn scalar(&self, state: FockState) -> i32 {
        // Scalar product initialization
        let mut scalar: i32 = 0;

        // Verifiying if orthogonal states
        if self.is_null || state.is_null {
            scalar = 0;
        } else if self.integer == state.integer {
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
    fn create(&mut self, index: u32) -> FockState {
        // Initializing new state binary number or status (null or not)
        let mut to_null: bool = false;
        let mut new_state: i32 = self.integer;

        // Verifiying if indexing isn't too long for Hilbert space
        if index >= (self.n_sites as u32).pow(2) {
            std::println!("Cannot create fermion at index: {}", index);
            std::process::exit(1);

        // Verifying if a fermion if already at position 'index' in state
        } else if self.integer_to_binary()[index as usize] == 1 || self.is_null {
            to_null = true;
            new_state = 0;

        // Updating Fock state integer after creating fermion
        } else {
            new_state += (2 as i32).pow(index);
        }
        FockState {
            n_sites: self.n_sites,
            integer: new_state,
            is_null: to_null,
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
    fn destroy(&mut self, index: u32) -> FockState {
        // Initializing new state binary number or status (null or not)
        let mut to_null: bool = false;
        let mut new_state: i32 = self.integer;

        // Verifiying if indexing isn't too long for Hilbert space
        if index >= (self.n_sites as u32).pow(2) {
            std::println!("Cannot create fermion at index: {}", index);
            std::process::exit(1);

        // Verifying if no fermions are at position 'index' in state
        } else if self.integer_to_binary()[index as usize] == 0 || self.is_null {
            to_null = true;
            new_state = 0;

        // Updating Fock state integer after destroying fermion
        } else {
            new_state -= (2 as i32).pow(index);
        }
        FockState {
            n_sites: self.n_sites,
            integer: new_state,
            is_null: to_null,
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
    fn number(&mut self, index: u32) -> FockState {
        // Initializing new state binary number or status (null or not)
        let mut to_null: bool = false;
        let mut new_state: i32 = self.integer;

        // Verifiying if indexing isn't too long for Hilbert space
        if index >= (self.n_sites as u32).pow(2) {
            std::println!("Fermion cannot be at index: {}", index);
            std::process::exit(1);

        // Verifiying if a fermion at site 'index' or if state is null
        } else if self.integer_to_binary()[index as usize] == 0 || self.is_null {
            new_state = 0;
            to_null = true;
        }
        FockState {
            n_sites: self.n_sites,
            integer: new_state,
            is_null: to_null,
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
    pub fn interaction_term(&self, state0: i32, state1: i32) -> i32 {
        // Initializing matrix element
        let mut coefficient: i32 = 0;

        // Initializing 'bra'
        let state0: FockState = FockState {
            n_sites: self.n_sites,
            integer: state0,
            is_null: false,
        };

        // Initializing 'ket'
        let mut state1: FockState = FockState {
            n_sites: self.n_sites,
            integer: state1,
            is_null: false,
        };

        // Main loop over number of sites in the cluster (i)
        for site in 0..self.n_sites as u32 {
            let mut inter_1: FockState = state1.number(site);
            let inter_2: FockState = inter_1.number(site + self.n_sites as u32);
            coefficient += self.u * state0.scalar(inter_2);
        }
        coefficient
    }

    pub fn kinetic_term(&self, state0: i32, state1: i32) -> i32 {
        // Initializing matrix element
        let mut coefficient: i32 = 0;

        // Initializing 'bra'
        let state0: FockState = FockState {
            n_sites: self.n_sites,
            integer: state0,
            is_null: false,
        };

        // Initializing 'ket'
        let mut state1: FockState = FockState {
            n_sites: self.n_sites,
            integer: state1,
            is_null: false,
        };

        // Main loop over number of sites in the cluster (i, j)
        for site_i in 0..self.n_sites as u32 {
            for site_j in 0..self.n_sites as u32 {
                // Removing 'on site' hoppings
                if site_i != site_j {
                    // Kinetic term for spin 'up'
                    let mut kin_1_up: FockState = state1.destroy(site_j);
                    let kin_2_up: FockState = kin_1_up.create(site_i);
                    coefficient += self.t * state0.scalar(kin_2_up);

                    // Kinetic term for spin 'down'
                    let mut kin_1_down: FockState = state1.destroy(site_j + self.n_sites as u32);
                    let kin_2_down: FockState = kin_1_down.create(site_i + self.n_sites as u32);
                    coefficient += self.t * state0.scalar(kin_2_down);
                }
            }
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
