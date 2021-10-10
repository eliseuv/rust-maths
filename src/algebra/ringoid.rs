//! Ringoids are ring-like algebraic structures.
//! They are algebraic structures built from a set closed under two binary operations, often called addition and multiplication, with multiplication being distributive over addition.
//!
//! Ringoid algebraic structures:
//! - Ringoid: set with two binary operations, often called addition and multiplication, with multiplication being distributive over addition
//! - Semiring: ringoid such that is also a monoid under each binary operation (the multiplicative identity is usually called 1 and the additive identity is usually called 0)
//! - Near ring: semiring whose addition monoid is a group
//! - Ring: semiring whose addition monoid is an abelain group
//! - Commutative ring: ring whose multiplication operation is commutative
//! - Field: commutative ring which contains a multiplicative inverse for every nonzero element (i.e., except for the identity of the addition monoid)

use super::BaseSet;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg};

//=========//
// Ringoid //
//=========//

// Ringoid is a set with two binary operations, often called addition and multiplication, with multiplication being distributive over addition
// ∀ a,b ∈ S, a∙b ∈ S, a+b ∈ S
// ∀ a,b,c ∈ S, a∙(b+c) = a∙b+a∙c
pub trait AbstractRingoid<Set> {
    fn add(lhs: Set, rhs: Set) -> Set;
    fn mult(lhs: Set, rhs: Set) -> Set;
}

// Ringoid trait object
pub trait Ringoid: Sized + Mul<Output = Self> + MulAssign + Add<Output = Self> + AddAssign {}

// Ringoid tait objects implement ringoid abstract trait
impl<Set: Ringoid> AbstractRingoid<Set> for Set {
    fn add(lhs: Set, rhs: Set) -> Set {
        lhs + rhs
    }

    fn mult(lhs: Set, rhs: Set) -> Set {
        lhs * rhs
    }
}

// The ringoid trait is automatically satisfied for every sized object that has the usual multiplication and addition operators defined
impl<Set> Ringoid for Set where
    Set: Sized + Mul<Output = Self> + MulAssign + Add<Output = Self> + AddAssign
{
}

//==========//
// Semiring //
//==========//

// Semiring is a ringoid such that S is a monoid under each of the two operations
// This implies the existence of an identity for each operation
// ∃! 1 ∈ S, ∀ a ∈ S, a∙1 = 1∙a = a
// ∃! 0 ∈ S, ∀ a ∈ S, a+0 = 0+a = a
pub trait AbstractSemiring<Set>: AbstractRingoid<Set> {
    // Multiplicative identity
    fn id_mult() -> Set;
    // Additive identity
    fn id_add() -> Set;
}

// Semiring trait object
pub trait Semiring: Ringoid + AbstractSemiring<Self> {
    // Get multiplicative identity
    fn get_id_mult() -> Self;
    // Set to multiplicative identity
    fn set_id_mult(&mut self) {
        *self = Self::id_mult();
    }
    // Test if is multiplicative identity
    fn is_id_mult(&self) -> bool;
    // Get additive identity
    fn get_id_add() -> Self;
    // Set to additive identity
    fn set_id_add(&mut self) {
        *self = Self::id_add();
    }
    // Test if is additive identity
    fn is_id_add(&self) -> bool;
}

// Semiring trait objects implement semiring abstract trait
impl<Set: Semiring> AbstractSemiring<Set> for Set {
    fn id_mult() -> Set {
        Self::get_id_mult()
    }
    fn id_add() -> Set {
        Self::get_id_add()
    }
}

//===========//
// Near Ring //
//===========//

// Near-ring is a semiring whose additive monoid is a group
// This implies that every element has an inverse with respect to addition
// ∀ a ∈ S, ∃! -a, a+(-a)=(-a)+a = 0
pub trait AbstractNearRing<Set>: AbstractSemiring<Set> {
    // Get additive inverse of a given element
    fn neg(a: Set) -> Set;
}

// Near ring trait object
pub trait NearRing: Ringoid + AbstractSemiring<Self> + Neg<Output = Self> {
    fn as_neg(self) -> Self {
        -self
    }
}

// Near ring trait object implements near ring abstract trait
impl<Set: NearRing> AbstractNearRing<Set> for Set {
    fn neg(a: Set) -> Set {
        a.as_neg()
    }
}

//======//
// Ring //
//======//

// Ring is a semiring whose additive monoid is an abelian group
// This implies that every element has an inverse with respect to addition and that addition is commutative
// ∀ a,b ∈ S, a+b=b+a
#[marker]
pub trait AbstractRing<Set>: AbstractNearRing<Set> {}

// Ring trait object
#[marker]
pub trait Ring: NearRing {}

// Ring trait object satisfies ring abstract trait
impl<Set: Ring> AbstractRing<Set> for Set {}

//==================//
// Commutative Ring //
//==================//

// Commutative ring is a ring whose multiplication operation is commutative
// ∀ a,b ∈ S, a·b=b·a
#[marker]
pub trait AbstractCommutativeRing<Set>: AbstractRing<Set> {}

// Commutative ring trait object
#[marker]
pub trait CommutativeRing: Ring {}

// Commutative ring trait object implements commutative ring abstract trait
impl<Set: CommutativeRing> AbstractCommutativeRing<Set> for Set {}

// Implement ring for commutative ring
impl<T, Set> AbstractRing<Set> for T where T: AbstractCommutativeRing<Set> {}
impl<Set: CommutativeRing> Ring for Set {}

//=======//
// Field //
//=======//

// Field is a commutative ring which contains a multiplicative inverse for every nonzero element
// ∀ a ∈ S, a≠0 ∈ S, ∃! a', a·a'=a'·a=1
pub trait AbstractField<Set>: AbstractCommutativeRing<Set> {
    // Get inverse
    fn inv(a: &Set) -> Result<Set, &str>;
}

// Field trait object
pub trait Field: CommutativeRing {
    fn inv(&self) -> Result<Self, &str>
    where
        Self: Sized;
}

// Field trait object implements field abstract trait
impl<Set: Field> AbstractField<Set> for Set {
    fn inv(a: &Set) -> Result<Set, &str> {
        a.inv()
    }
}

// Implement commutative ring for field
impl<T, Set> AbstractCommutativeRing<Set> for T where T: AbstractField<Set> {}
impl<Set: Field> CommutativeRing for Set {}
