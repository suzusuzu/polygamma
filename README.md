# polygamma

![Rust](https://github.com/suzusuzu/polygamma/workflows/Rust/badge.svg) 
[![Document](https://docs.rs/polygamma/badge.svg)](https://docs.rs/polygamma) 
[![Crates](https://img.shields.io/crates/v/polygamma.svg)](https://crates.io/crates/polygamma)

This library provides Polygamma function.
It is forked from [http://www.rd.dnc.ac.jp/~tunenori/src/polygamma.c](http://www.rd.dnc.ac.jp/~tunenori/src/polygamma.c)

# Usage

```rust
use polygamma::polygamma;

fn main() {
    println!("digamma(1.0)  = {}", polygamma(0, 1.0).unwrap());
    println!("trigamma(1.0) = {}", polygamma(1, 1.0).unwrap());
}
```
