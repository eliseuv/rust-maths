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

// TODO: Use traits defined in the groupoid module.

//=========//
// Ringoid //
//=========//

// Here we use the same name for the property and the simplest structure that implements it.
// Ringoid: set with two binary operations, often called <<addition>> and <<multiplication>>, with <<multiplication>> being distributive over <<addition>>.
// ∀ a,b ∈ S, a∙b ∈ S, a+b ∈ S
// ∀ a,b,c ∈ S, a∙(b+c) = a∙b+a∙c

pub trait AbstractRingoid<Set: BaseSet> {
    fn add(lhs: &Set, rhs: &Set) -> Set;
    fn mult(lhs: &Set, rhs: &Set) -> Set;
}

pub trait Ringoid: BaseSet {
    fn add(&self, other: &Self) -> Self;
    fn mult(&self, other: &Self) -> Self;
}

impl<Set: Ringoid> AbstractRingoid<Set> for Set {
    fn add(lhs: &Set, rhs: &Set) -> Set {
        lhs.add(rhs)
    }
    fn mult(lhs: &Set, rhs: &Set) -> Set {
        lhs.mult(rhs)
    }
}

//=============================//
// Commutativity of Operations //
//=============================//

// Commutativity of <<addition>>
// ∀ a,b ∈ S, a+b = b+a
#[marker]
pub trait AbstractRingoidCommutativeAddition<Set: BaseSet>: AbstractRingoid<Set> {}
#[marker]
pub trait RingoidCommutativeAddition: Ringoid {}
impl<Set: RingoidCommutativeAddition> AbstractRingoidCommutativeAddition<Set> for Set {}

// Commutativity of <<Multiplication>>
// ∀ a,b ∈ S, a·b = b·a
#[marker]
pub trait AbstractRingoidCommutativeMultiplication<Set: BaseSet>: AbstractRingoid<Set> {}
#[marker]
pub trait RingoidCommutativeMultiplication: Ringoid {}
impl<Set: RingoidCommutativeMultiplication> AbstractRingoidCommutativeMultiplication<Set> for Set {}

//====================//
// Identity Elemenets //
//====================//

// Existence of identity elements for each of the two ringoid's binary operations.
// Often called 1 for the <<multiplication>> and 0 for the <<addition>>.
// ∃! 1 ∈ S, ∀ a ∈ S, a∙1 = 1∙a = a
// ∃! 0 ∈ S, ∀ a ∈ S, a+0 = 0+a = a

pub trait AbstractRingoidIdentities<Set: BaseSet>: AbstractRingoid<Set> {
    fn id_add() -> Set;
    fn id_mult() -> Set;
}

pub trait RingoidIdentities: Ringoid {
    fn id_add() -> Self;
    fn set_id_add(&mut self) {
        *self = Self::id_add();
    }
    fn is_id_add(&self) -> bool {
        self == &Self::id_add()
    }
    fn id_mult() -> Self;
    fn set_id_mult(&mut self) {
        *self = Self::id_mult();
    }
    fn is_id_mult(&self) -> bool {
        self == &Self::id_add()
    }
}

impl<Set: RingoidIdentities> AbstractRingoidIdentities<Set> for Set {
    fn id_add() -> Set {
        Set::id_add()
    }
    fn id_mult() -> Set {
        Set::id_mult()
    }
}

//==================================//
// Inverse Element for <<Addition>> //
//==================================//

// Inverse of any element with respect to <<addition>>.
// ∀ a ∈ S, ∃! -a, (a)+(-a)=(-a)+(a) = 0

pub trait AbstractRingoidNegative<Set: BaseSet>: AbstractRingoidIdentities<Set> {
    fn neg(a: &Set) -> Set;
}

pub trait RingoidNegative: RingoidIdentities {
    fn neg(&self) -> Self;
}

impl<Set: RingoidNegative> AbstractRingoidNegative<Set> for Set {
    fn neg(a: &Set) -> Set {
        a.neg()
    }
}

//========================================//
// Inverse Element for <<Multiplication>> //
//========================================//

// Inverse of any NONZERO element with respect to <<multiplication>>.
// ∀ a ∈ S, a≠0 ∈ S, ∃! a', a·a'=a'·a=1

pub trait AbstractRingoidInverse<Set: BaseSet>: AbstractRingoidIdentities<Set> {
    fn inv(a: &Set) -> Result<Set, &str>;
}

pub trait RingoidInverse: RingoidIdentities {
    fn inv(&self) -> Result<Self, &str>;
}

impl<Set: RingoidInverse> AbstractRingoidInverse<Set> for Set {
    fn inv(a: &Set) -> Result<Set, &str> {
        a.inv()
    }
}

//============//
// Structures //
//============//

//==========//
// Semiring //
//==========//

// Semiring: ringoid such that S is a monoid under each of the two operations

pub trait AbstractSemiring<Set: BaseSet> = AbstractRingoid<Set> + AbstractRingoidIdentities<Set>;

pub trait Semiring = Ringoid + RingoidIdentities;

//===========//
// Near-Ring //
//===========//

// Near-ring: semiring whose additive monoid is a group

pub trait AbstractNearRing<Set: BaseSet> = AbstractSemiring<Set> + AbstractRingoidNegative<Set>;

pub trait NearRing = Semiring + RingoidNegative;

//======//
// Ring //
//======//

// Ring: semiring whose additive monoid is an abelian group

pub trait AbstractRing<Set: BaseSet> =
    AbstractNearRing<Set> + AbstractRingoidCommutativeAddition<Set>;

pub trait Ring = NearRing + RingoidCommutativeAddition;

//==================//
// Commutative Ring //
//==================//

// Commutative ring: ring whose multiplication operation is commutative

pub trait AbstractCommutativeRing<Set: BaseSet> =
    AbstractRing<Set> + AbstractRingoidCommutativeMultiplication<Set>;

pub trait CommutativeRing = Ring + RingoidCommutativeMultiplication;

//=======//
// Field //
//=======//

// Field: commutative ring which contains a multiplicative inverse for every nonzero element

pub trait AbstractField<Set: BaseSet> = AbstractCommutativeRing<Set> + AbstractRingoidInverse<Set>;

pub trait Field = CommutativeRing + RingoidInverse;
