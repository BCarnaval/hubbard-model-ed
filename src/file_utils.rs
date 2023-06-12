// Contains useful functions to generate and manage csv data files.
// This module is using the 'csv' crate.

use indicatif::{ProgressBar, ProgressStyle};
use std::error::Error;
use std::fs::OpenOptions;

/// Initialize a ProgressBar using 'indicatif' crate.
///
/// # Examples
///
/// ```rust
/// let msg: String = String::from("A message");
/// let iters: u64 = 1000;
/// let pb = init_progress_bar(&msg, iters);
/// ```
pub fn init_progress_bar(msg: String, iters: u64) -> ProgressBar {
    let pb = ProgressBar::new(iters);
    pb.set_style(
        ProgressStyle::with_template(&format!(
            "[{}] {{elapsed_precise}} {{bar:40.cyan/blue}} {{pos:>2}}/{{len:2}}",
            &msg
        ))
        .unwrap()
        .progress_chars("=>-"),
    );
    pb
}

/// Returns text file number of row(s).
///
/// # Examples
///
/// ```rust
/// let path = String::from("./path/to/model");
/// let nrows = file_lenght(&path).unwrap();
/// println!("File {} has {} lines", &path, nrows);
/// ```
#[allow(dead_code)]
fn file_length(file_path: &String, headers: bool) -> Result<u64, Box<dyn Error>> {
    let rdr = csv::ReaderBuilder::new()
        .has_headers(headers)
        .from_path(file_path)
        .unwrap();

    let mut recs = rdr.into_records();
    let mut nrows: u64 = 0;
    loop {
        let next_row = recs.reader().position().line();
        if recs.next().is_none() {
            break;
        }
        nrows = next_row;
    }
    Ok(nrows) // Removing headers row
}

/// Initialize file writter using 'csv' crate.
///
/// # Examples
///
/// ```rust
/// let path = String::from("path/to/file");
/// let writter = init_file_writter(&path, true);
/// ```
pub fn init_file_writter(path: &String, has_headers: bool) -> csv::Writer<std::fs::File> {
    let file_writter = OpenOptions::new()
        .truncate(true)
        .write(true)
        .create(true)
        .open(&path)
        .unwrap();

    let writter = csv::WriterBuilder::new()
        .flexible(true)
        .has_headers(has_headers)
        .delimiter(b' ')
        .from_writer(file_writter);

    writter
}
