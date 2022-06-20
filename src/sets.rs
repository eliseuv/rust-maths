//! Basic set theory concepts

/// Base set
///
/// A set must have the size of its elements known at compile time and have elements distinguishable from one another.
pub trait BaseSet = Sized + PartialEq;

/// Finite set
///
/// A finite set is a base set with a finite known number of elements.
pub trait FiniteSet: BaseSet {
    /// Get the number of elements in this set
    fn element_count(&self) -> usize;
}
