<div align="center">

# hubbard_model_ed
  
This repository contains an exact diagonalization implementation for the Hubbard model (in the approximation of 1D spins chain with first neigbohrs hopping terms using periodic boundary conditions) defined as
$$H = H_t + H_U = \sum_{i, j, \sigma}t_{ij}(c^\dagger_{i\sigma}c_{j\sigma} + h.c.) + U\sum_i n_{i\uparrow}n_{i\downarrow},$$
where $c^\dagger$ and $c$ are respectively the second quantization creation/anihilation operators and where $n$ represents the number operator from the same formalism. The code is entirely written in [Rust](https://www.rust-lang.org/) and will eventually be parallelized using [rayon](https://github.com/rayon-rs/rayon) Rust crate.

![Rust](https://img.shields.io/badge/rust-%23000000.svg?style=for-the-badge&logo=rust&logoColor=white)
  
![Build](https://img.shields.io/github/actions/workflow/status/BCarnaval/hubbard_model_ed/rust.yml?color=%23a3d1af&style=for-the-badge) ![LICENSE](https://img.shields.io/github/license/BCarnaval/hubbard_model_ed?color=blue&style=for-the-badge) ![release](https://img.shields.io/github/v/tag/BCarnaval/hubbard_model_ed?color=%23FF7F50&style=for-the-badge)
  
</div>

## Table of contents

[Requirements](#requirements)

[Installation](#installation)

[Usage](#usage)

[TODO](#todo)

## Requirements

## Installation

## Usage

## TODO

- [ ] Complete the `README.md`
- [x] Include periodic boundary conditions to hoppings operator (phase correction)
- [ ] Include parallel computing using [rayon](https://github.com/rayon-rs/rayon)
- [x] Comment the code base
- [x] Unit testing
