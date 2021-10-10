//! Gropoids are group-like algebraic structures.
//! They are algebraic structures built from a set closed under a binary operation.
//!
//! Groupoid algebraic structures:
//! - Magma: set closed under a binary operation
//! - Semigroup: associative magma
//! - Monoid: semigroup with unit element
//! - Group: monoid with an inverse for every element
//! - Abelian group: commutative group

use super::BaseSet;
use std::ops::{Mul, MulAssign};

//==================//
// Binary Operation //
//==================//

// Binary operation defined over a given set
pub trait BinaryOperation<Set: BaseSet> {
    // op: S ⨉ S → S
    //     op(a,b) ↦ a·b
    // ∀ a,b ∈ S, a∙b ∈ S
    fn op(lhs: Set, rhs: Set) -> Set;
}

//==============================//
// Associative Binary Operation //
//==============================//

// Associative binary operation defined over a given set
// ∀ a,b,c ∈ S, (a∙b)∙c = a∙(b∙c)
#[marker]
pub trait AssociativeBinaryOperation<Set: BaseSet>: BinaryOperation<Set> {}

//==============================//
// Commutative Binary Operation //
//==============================//

// Commutative binary operation defined over a given set
// ∀ a,b ∈ S, a∙b = b∙a
#[marker]
pub trait CommutativeBinaryOperation<Set: BaseSet>: BinaryOperation<Set> {}

// TODO: Auto implement BinOp if CommBinOp is implemented without collision
// If associative binary operation trait is satified so is binary operation
/*
impl<T, Set> BinaryOperation<Set> for T
where
    Set: BaseSet,
    T: AssociativeBinaryOperation<Set>,
{
    fn op(lhs: Set, rhs: Set) -> Set {
        AssociativeBinaryOperation::<Set>::op(lhs, rhs)
    }
}
 */

//==================//
// Identity Element //
//==================//

// Defines the identity element of a given set
pub trait IdentityElement<Set: BaseSet>: BinaryOperation<Set> {
    // ∃! e ∈ S, ∀ a ∈ S, a∙e = e∙a = a
    fn id() -> Set;
}

//=================//
// Inverse Element //
//=================//

// Returns the inverse of ANY given element in a set
pub trait InverseElement<Set: BaseSet>: BinaryOperation<Set> + IdentityElement<Set> {
    // ∀ a ∈ S, ∃! a' ∈ S, a∙a' = a'∙a = e
    fn inv(a: Set) -> Set;
}

//=======//
// Magma //
//=======//

// Magma: set closed under a binary operation

// Abstract magma structure over a given set
#[marker]
pub trait AbstractMagma<Set: BaseSet>: BinaryOperation<Set> {}

// All types defining the required properties, satisfy abstract trait
impl<T, Set> AbstractMagma<Set> for T
where
    Set: BaseSet,
    T: BinaryOperation<Set>,
{
}

// Magma trait object is a set that implements a binary operation
pub trait Magma: BaseSet + BinaryOperation<Self> {
    fn op(self, other: Self) -> Self;
}

// Once trait object is implemented, also implement abstract properties
impl<Set: Magma> BinaryOperation<Set> for Set {
    fn op(lhs: Set, rhs: Set) -> Set {
        lhs.op(rhs)
    }
}

// Trait object implements abstract trait
impl<Set: Magma> AbstractMagma<Set> for Set {}

// Implement magma trait for sets with multiplication defined
impl<Set> Magma for Set
where
    Set: BaseSet + Mul<Output = Self> + MulAssign,
{
    fn op(self, other: Self) -> Self {
        self * other
    }
}

//===========//
// Semigroup //
//===========//

// Semigroup: associative magma

// Abstract semigroup over a given set
#[marker]
pub trait AbstractSemigroup<Set: BaseSet>: AssociativeBinaryOperation<Set> {}

// All types defining the required properties, satisfy abstract trait
impl<T, Set> AbstractSemigroup<Set> for T
where
    Set: BaseSet,
    T: AssociativeBinaryOperation<Set>,
{
}

// Semigroup trait object is a set that implements an associative binary operation
#[marker]
pub trait Semigroup: Magma + AssociativeBinaryOperation<Self> {}

// Once trait object is implemented, also implement abstract properties
impl<Set: Semigroup> AssociativeBinaryOperation<Set> for Set {}

// Trait object implements abstract trait
impl<Set: BaseSet + Semigroup> AbstractSemigroup<Set> for Set {}

//========//
// Monoid //
//========//

// Monoid: semigroup with identity element

// Abstract monoid over a given set
#[marker]
pub trait AbstractMonoid<Set: BaseSet>: AbstractSemigroup<Set> + IdentityElement<Set> {}

// All types defining the required properties, satisfy abstract trait
impl<T, Set> AbstractMonoid<Set> for T
where
    Set: BaseSet,
    T: AbstractSemigroup<Set> + IdentityElement<Set>,
{
}

// Monoid trait object
pub trait Monoid: Semigroup + IdentityElement<Self> {
    // Get identity
    fn id() -> Self;
    // Set to identity
    fn set_id(&mut self) {
        *self = Monoid::id();
    }
    // Test if is identity
    fn is_id(&self) -> bool;
}

// Once trait object is implemented, also implement abstract properties
impl<Set: Monoid> IdentityElement<Set> for Set {
    fn id() -> Set {
        Monoid::id()
    }
}

// Trait object implements abstract trait
impl<Set: Monoid> AbstractMonoid<Set> for Set {}

//====================//
// Commutative monoid //
//====================//

// Commutative monoid: monoid whose binary operation is also commutative

// Abstract commutative monoid over a set
#[marker]
pub trait AbstractCommutativeMonoid<Set: BaseSet>:
    AbstractMonoid<Set> + CommutativeBinaryOperation<Set>
{
}

// All types defining the required properties, satisfy abstract trait
impl<T, Set> AbstractCommutativeMonoid<Set> for T
where
    Set: BaseSet,
    T: AbstractMonoid<Set> + CommutativeBinaryOperation<Set>,
{
}

// Commutative monoid trait object
#[marker]
pub trait CommutativeMonoid: Monoid + CommutativeBinaryOperation<Self> {}

// Once trait object is implemented, also implement abstract properties
impl<Set: CommutativeMonoid> CommutativeBinaryOperation<Set> for Set {}

// Trait object implements abstract trait
impl<Set: CommutativeMonoid> AbstractCommutativeMonoid<Set> for Set {}

//=======//
// Group //
//=======//

// Group: monoid for which every element has an inverse

// Abstract group over a set
#[marker]
pub trait AbstractGroup<Set: BaseSet>: AbstractMonoid<Set> + InverseElement<Set> {}

// All types defining the required properties, satisfy abstract trait
impl<T, Set> AbstractGroup<Set> for T
where
    Set: BaseSet,
    T: AbstractMonoid<Set> + InverseElement<Set>,
{
}

// Group trait object
pub trait Group: Monoid + InverseElement<Self> {
    fn inv(&self) -> Self;
}

// Once trait object is implemented, also implement abstract properties
impl<Set: Group> InverseElement<Set> for Set {
    fn inv(a: Set) -> Set {
        a.inv()
    }
}

// Trait object implements abstract trait
impl<Set: Group> AbstractGroup<Set> for Set {}

//===============//
// Abelian Group //
//===============//

// Abelian group: group whose binary operation is commutative

// Abstract abelian group over a set
#[marker]
pub trait AbstractAbelianGroup<Set: BaseSet>:
    AbstractCommutativeMonoid<Set> + InverseElement<Set>
{
}

// All types defining the required properties, satisfy abstract trait
impl<Set, T> AbstractAbelianGroup<Set> for T
where
    Set: BaseSet,
    T: AbstractCommutativeMonoid<Set> + InverseElement<Set>,
{
}

// Abelian group trait object
#[marker]
pub trait AbelianGroup: CommutativeMonoid + InverseElement<Self> {}

// Abelian group trait object also implements abstract trait
impl<Set: AbelianGroup> AbstractAbelianGroup<Set> for Set {}

pub mod instances {
    //! Example implementations of groupoids
    //!

    use super::{BinaryOperation, IdentityElement, InverseElement};

    // Cyclic group of order 2 (Z_2)
    // U ∙ A = A ∙ U = A
    // A ∙ A = U
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub enum Z2 {
        U,
        A,
    }

    impl BinaryOperation<Self> for Z2 {
        fn op(lhs: Self, rhs: Self) -> Self {
            match lhs {
                Z2::U => match rhs {
                    Z2::U => Z2::U,
                    Z2::A => Z2::A,
                },
                Z2::A => match rhs {
                    Z2::U => Z2::A,
                    Z2::A => Z2::U,
                },
            }
        }
    }

    impl IdentityElement<Self> for Z2 {
        fn id() -> Self {
            Z2::U
        }
    }

    impl InverseElement<Self> for Z2 {
        fn inv(a: Self) -> Self {
            match a {
                Z2::U => Z2::U,
                Z2::A => Z2::A,
            }
        }
    }
}

// Tests
#[cfg(test)]
mod tests {

    use super::instances::*;

    #[test]
    fn groups() {
        let u = Z2::U;
        let a = Z2::A;
        // Test binary operation
        // Test unity
        // Test inverses
    }
}
