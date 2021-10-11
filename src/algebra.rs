//! Simple representation of ususal algebraic structures using Rust's trait system
//!
//! `AbstractStrucuture<Set>` are marker traits that indicates that a given type has a given algebraic structure's properties defined over a given set (type).
//! `Structure` are trait objects equiped with a given algebraic structure. Here we call the set (type) itself the algebraic structure.
//!
//! Every trait object also implements its corresponding abstract trait.
//!
//! Since the algebraic structure is completely defined by its properties, the structure trait is defined as type alias for its properties.
//! Once the properties are implemented the type is marked as given algebric structure.

//==========//
// Base Set //
//==========//

// A set must have the size of its elements known at compile time and have elements distinguishable from one another
pub trait BaseSet = Sized + PartialEq;

// Group-like structures
pub mod groupoid;
// Ring-like structures
pub mod ringoid;
