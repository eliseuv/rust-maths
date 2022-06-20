//! Physica is a scientific library written as personal project for learning the Rust language.

#![feature(marker_trait_attr)]
#![feature(trait_alias)]

// Set theory concepts
pub mod sets;
// Abstract algebra
pub mod algebra;
// Topology
pub mod topology;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
