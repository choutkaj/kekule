use std::{error::Error, process};

mod compare;
mod dataset;
mod features;
mod observation;
mod runner;

fn boxed_error(message: impl Into<String>) -> Box<dyn Error> {
    std::io::Error::other(message.into()).into()
}

fn main() {
    if let Err(error) = runner::run() {
        eprintln!("error: {error}");
        process::exit(1);
    }
}

#[cfg(test)]
mod tests;
