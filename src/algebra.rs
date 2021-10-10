//! Simple representation of ususal algebraic structures using Rust's trait system
//!
//! `AbstractStrucuture<Set>` are marker traits that indicates that a given type has a given algebraic structure's properties defined over a given set (type).
//! `Structure` are trait objects equiped with a given algebraic structure. Here we call the set (type) itself the algebraic structure.
//!
//! The algebraic structure is completely defined by its properties.
//! The abstract trait requires the implemenentation of the properties, but also the implementation of such properties implies that trait is satisfied.
//! Once the properties are implemented the type is marked as given algebric structure.
//!
//! Every trait object also implements its corresponding abstract trait.

//==========//
// Base Set //
//==========//

// A set must have the size of its elements known at compile time and have elements distinguishable from one another
#[marker]
pub trait BaseSet: Sized + PartialEq {}

// Every type that satisfies these conditions can be considered a set
impl<Set> BaseSet for Set where Set: Sized + PartialEq {}

// Group-like structures
pub mod groupoid;
// Ring-like structures
pub mod ringoid;
