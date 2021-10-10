//! Physica is a scientific library written as personal project for learning the Rust language.

#![feature(marker_trait_attr)]
//#![feature(min_specialization)]

// Abstract algebra
pub mod algebra;
// Topology
pub mod topology {}
// Linear algebra
pub mod linalg {}
// Matrices numerical methods
pub mod matrix {}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
