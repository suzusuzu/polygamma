//!
//! This library provides Polygamma function.
//! It is forked from [http://www.rd.dnc.ac.jp/~tunenori/src/polygamma.c](http://www.rd.dnc.ac.jp/~tunenori/src/polygamma.c)
//!
//! # Usage
//!
//! ```
//! use polygamma::polygamma;
//! fn main() {
//!     println!("digamma(1.0)  = {}", polygamma(0, 1.0).unwrap());
//!     println!("trigamma(1.0) = {}", polygamma(1, 1.0).unwrap());
//! }
//! ```
//!

mod polygamma_;
pub use polygamma_::polygamma;