// This is a comment, and is ignored by the compiler.
// You can test this code by clicking the "Run" button over there ->
// or if you prefer to use your keyboard, you can use the "Ctrl + Enter"
// shortcut.

use itertools::Itertools;
use std::{println, vec};

use crate::{array_utils::lapack_diagonalization, file_utils::init_file_writter};

#[derive(Debug)]
pub struct FockState {
    // Public attributes
    pub n_sites: u32,
    pub integer: i32,

    // Private attributes
    is_null: bool,
    sign: i32,
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
    fn scalar(&self, state: i32) -> f32 {
        // Scalar product initialization
        let mut scalar: f32 = 0.;

        // Verifiying if orthogonal states
        if self.is_null {
            scalar = 0.;
        } else if self.integer == state {
            scalar = 1.;
        }
        scalar
    }

    /// Second quantization creation operator definition.
    ///
    /// Examples
    ///
    /// ```rust
    /// let state0 = FockState { n_sites: 2, integer: 5, is_null: false , sign: 1 };
    /// state0.create(1);
    /// println!("{:?}", state0.integer_to_binary());
    /// ```
    fn create(&mut self, index: i32) {
        // Defining binary mask
        let mask: i32 = 1 << ((2 * self.n_sites as i32 - 1) - index);
        let mut abs_state: i32 = self.integer.abs();

        // Verifying if a fermion if already at position 'index' in state
        if abs_state & mask != 0 || self.is_null {
            self.is_null = true;

        // Updating Fock state integer after creating fermion
        } else {
            abs_state ^= mask;
        }
        self.integer = self.sign * abs_state
    }

    /// Second quantization anihilation operator definition.
    ///
    /// Examples
    ///
    /// ```rust
    /// let state0 = FockState { n_sites: 2, integer: 5, is_null: false , sign: 1 };
    /// state0.destroy(1);
    /// println!("{:?}", state0.integer_to_binary());
    /// ```
    fn destroy(&mut self, index: i32) {
        // Defining binary mask
        let mask: i32 = 1 << ((2 * self.n_sites as i32 - 1) - index);
        let mut abs_state: i32 = self.integer.abs();

        // Verifying if no fermions are at position 'index' in state
        if abs_state & mask == 0 || self.is_null {
            self.is_null = true;

        // Updating Fock state integer after destroying fermion
        } else {
            abs_state ^= mask;
        }
        self.integer = self.sign * abs_state
    }

    /// Second quantization number operator definition.
    ///
    /// Examples
    ///
    /// ```rust
    /// let state0 = FockState { n_sites: 2, integer: 5, is_null: false , sign: 1 };
    /// state0.number(1);
    /// println!("{:?}", state0.integer_to_binary());
    /// ```
    fn number(&mut self, index: u32) {
        // Verifiying if a fermion at site 'index' or if state is null
        let mask: i32 = 1 << ((2 * self.n_sites - 1) - index);
        if self.integer & mask == 0 || self.is_null {
            self.is_null = true;
        }
    }
}

#[derive(Debug)]
pub struct Hubbard {
    // Public attributes
    pub n_sites: u32,
    pub t: f32,
    pub u: f32,
}

impl Hubbard {
    /// Computes 'on-site' interaction for given Fock State using second
    /// quantization number operator.
    pub fn interaction_term(&self, state_0: i32) -> f32 {
        // Initializing matrix element
        let mut coefficient: f32 = 0.;

        // Main loop over number of sites in the cluster (i)
        for site in 0..self.n_sites {
            // Initializing 'ket'
            let mut ket_state: FockState = FockState {
                n_sites: self.n_sites,
                integer: state_0,
                is_null: false,
                sign: 1,
            };

            // Computing 'on site' interaction using number operator
            ket_state.number(site);
            ket_state.number(site + self.n_sites);

            // Updating matrix element value
            coefficient += self.u * ket_state.scalar(state_0);
        }
        coefficient
    }

    /// Computes first neighbours hoppings for given Fock State using second
    /// quantization operators.
    ///
    /// It outputs a vector containing linked states for given initial Fock state.
    pub fn kinetic_term(&self, state_0: i32) -> Vec<i32> {
        // Initializing subspace states
        let mut sub_states: Vec<i32> = Vec::new();
        let sites: Vec<i32> = (0..self.n_sites as i32).collect();

        // Main loop over number of sites in the cluster (i, j)
        for (idx, perm) in sites.iter().permutations(2).unique().enumerate() {
            // Defining neighbours coordinates using permutations
            let site_i: i32 = *perm[0];
            let site_j: i32 = *perm[1];
            let bound: i32 = self.n_sites as i32 - 1;

            // Initializing 'kets' (spin up & down)
            let mut ket_up: FockState = FockState {
                n_sites: self.n_sites,
                integer: state_0,
                is_null: false,
                sign: state_0.signum(),
            };

            let mut ket_down: FockState = FockState {
                n_sites: self.n_sites,
                integer: state_0,
                is_null: false,
                sign: state_0.signum(),
            };

            // Kinetic term (spin up & down)
            let ket_up_int: i32 = ket_up.integer;
            ket_up.destroy(site_j);
            ket_up.create(site_i);

            let ket_down_int: i32 = ket_down.integer;
            ket_down.destroy(site_j + self.n_sites as i32);
            ket_down.create(site_i + self.n_sites as i32);

            // Push new state inside subspace if not already there
            if !sub_states.contains(&ket_up.integer) && !ket_up.is_null {
                let mut sign_up: i32 = 1;
                // Verifying if current sites correspond to boundary conditions
                if idx == (bound - 1) as usize || idx == (bound * bound) as usize {
                    // Spin up phase adjustement
                    let ket_up_boundary: u32 = (ket_up_int >> self.n_sites).count_ones() - 1;
                    if !(ket_up_boundary % 2 == 0) {
                        sign_up = -1
                    }
                }
                sub_states.push(sign_up * ket_up.integer)
            }
            if !sub_states.contains(&ket_down.integer) && !ket_down.is_null {
                let mut sign_down: i32 = 1;
                // Verifying if current sites correspond to boundary conditions
                if idx == (bound - 1) as usize || idx == (bound * bound) as usize {
                    // Spin down phase ajustement
                    let down_temp_int: i32 = ket_down_int;
                    let shifted_int: i32 =
                        (down_temp_int >> self.n_sites as i32) << self.n_sites as i32;
                    let ket_down_boundary: u32 = (ket_down_int - shifted_int).count_ones() - 1;
                    if !(ket_down_boundary % 2 == 0) {
                        sign_down = -1
                    }
                }
                sub_states.push(sign_down * ket_down.integer)
            }
        }
        sub_states.sort_by_key(|i| i.abs());
        println!("{:?}", sub_states);
        sub_states
    }

    /// Finds a block of the Hubbard hamiltonian using one Fock State at a time.
    ///
    /// It outputs a tuple in which we found all the states involved in the current
    /// block  and the matrix element of the block sorted 'column-wise' as LAPACK would
    /// recommend.
    pub fn find_sub_block(&self, state: i32) -> (Vec<i32>, Vec<f32>) {
        // Test index for new substates
        let mut idx: u32 = 0;
        let mut sub_states: Vec<i32> = vec![state];

        // Matrix elements array
        let mut elems: Vec<f32> = Vec::new();

        // Continue loop until substates aren't new
        while idx < sub_states.len() as u32 {
            // Defining current sub state
            let current_state: i32 = sub_states[idx as usize];
            println!("{}", current_state);

            // Find first hopping states
            let new_states: Vec<i32> = self.kinetic_term(current_state);
            let mut filtered: Vec<i32> = new_states.clone();
            filtered.retain(|i: &i32| !sub_states.contains(&i.abs()));
            sub_states.append(&mut filtered);

            // Kinetic terms
            let states_copy: Vec<&i32> =
                sub_states.iter().filter(|i| i < &&current_state).collect();
            for sub_state in states_copy {
                if !new_states.contains(sub_state) {
                    elems.push(0.);
                } else {
                    elems.push(self.t);
                }
            }

            // On-site interaction coefficient
            elems.push(self.interaction_term(current_state));

            // Updating array parser
            sub_states.sort_by_key(|i| i.abs());
            idx += 1;
        }
        println!("\n{:?}", sub_states);
        (sub_states, elems)
    }

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
}

#[cfg(test)]
mod tests {

    use std::assert_eq;

    use crate::fock_space::{FockState, Hubbard};

    #[test]
    fn test_fock_scalar() {
        // Test Fock state has form: |1 1 ; 1 1>
        let test_state: FockState = FockState {
            n_sites: 2,
            integer: 15,
            is_null: false,
            sign: 1,
        };
        // Testing for null scalar
        let coefficient_0: f32 = test_state.scalar(0);
        assert_eq!(0., coefficient_0);

        // Testing for non zero scalar
        let coefficient_1: f32 = test_state.scalar(15);
        assert_eq!(1., coefficient_1);
    }

    #[test]
    fn test_fock_create() {
        // Test Fock state has form: |0 1 ; 1 1>
        let mut test_state: FockState = FockState {
            n_sites: 2,
            integer: 7,
            is_null: false,
            sign: 1,
        };
        // Testing for valid fermion creation
        // Answer should be: |1 1 ; 1 1>
        test_state.create(0);
        assert_eq!(15, test_state.integer);

        // Testing invalid fermion creation
        test_state.create(0);
        assert_eq!(true, test_state.is_null);
    }

    #[test]
    fn test_fock_destroy() {
        // Test Fock state has form: |0 1 ; 1 1>
        let mut test_state: FockState = FockState {
            n_sites: 2,
            integer: 7,
            is_null: false,
            sign: -1,
        };
        // Testing for valid fermion anihilation
        // Answer should be: |0 0 ; 1 1>
        test_state.destroy(1);
        assert_eq!(-3, test_state.integer);

        // Testing for invalid fermion anihilation
        test_state.destroy(1);
        assert_eq!(true, test_state.is_null);
    }

    #[test]
    fn test_fock_number() {
        // Test Fock state has form: |1 1 ; 0 1>
        let mut test_state: FockState = FockState {
            n_sites: 2,
            integer: 13,
            is_null: false,
            sign: 1,
        };
        // Testing for existing fermion
        test_state.number(0);
        assert_eq!(13, test_state.integer);

        // Testing for non-existing fermion
        test_state.destroy(0);
        test_state.number(0);
        assert_eq!(true, test_state.is_null);
    }

    #[test]
    fn test_hubbard_interaction() {
        // Test hubbard instance
        let test_model: Hubbard = Hubbard {
            n_sites: 2,
            t: 1.,
            u: 2.,
        };
        assert_eq!(4., test_model.interaction_term(15));
        assert_eq!(2., test_model.interaction_term(5));
        assert_eq!(0., test_model.interaction_term(1));
    }

    #[test]
    fn test_hubbard_hoppings() {
        // Test hubbard instance
        let test_model: Hubbard = Hubbard {
            n_sites: 2,
            t: 1.,
            u: 2.,
        };
        let empty: Vec<i32> = Vec::new();
        assert_eq!(empty, test_model.kinetic_term(0));
        assert_eq!(vec![2], test_model.kinetic_term(1));
        assert_eq!(vec![6, 9], test_model.kinetic_term(5));
    }

    #[test]
    fn test_hubbard_blocks() {
        // Test hubbard instance
        let test_model: Hubbard = Hubbard {
            n_sites: 3,
            t: 1.,
            u: 2.,
        };
        let sub_states: Vec<i32> = vec![9, 10, 12, 17, 18, 20, 33, 34, 36];
        let elements: Vec<f32> = vec![
            2., 1., 0., 1., 1., 0., 1., 0., 0., 0., 0., 1., 0., 1., 2., 0., 0., 1., 1., 1., 0., 1.,
            0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 1., 0., 1., 0., 0., 0., 1., 0., 0., 1., 1., 1.,
            2.,
        ];
        let (states, elems): (Vec<i32>, Vec<f32>) = test_model.find_sub_block(9);
        let empty: Vec<i32> = Vec::new();
        let difference: Vec<_> = sub_states
            .into_iter()
            .filter(|elem| !states.contains(&elem))
            .collect();
        assert_eq!(empty, difference);
        assert_eq!(elements, elems);
    }
}
