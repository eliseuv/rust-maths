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

//============//
// Properties //
//============//

//==================//
// Binary Operation //
//==================//

// A binary operation takes two members of a set and returns a member of the same set.
// op: S ⨉ S → S
//     op(a,b) ↦ a·b
// ∀ a,b ∈ S, a∙b ∈ S

pub trait AbstractBinaryOperation<Set: BaseSet> {
    fn op(lhs: &Set, rhs: &Set) -> Set;
}

pub trait BinaryOperation: BaseSet {
    fn op(&self, other: &Self) -> Self;
}

impl<Set: BinaryOperation> AbstractBinaryOperation<Set> for Set {
    fn op(lhs: &Set, rhs: &Set) -> Set {
        lhs.op(rhs)
    }
}

// Associativity
// ∀ a,b,c ∈ S, (a∙b)∙c = a∙(b∙c)

#[marker]
pub trait AbstractAssociativeBinaryOperation<Set: BaseSet>: AbstractBinaryOperation<Set> {}

#[marker]
pub trait AssociativeBinaryOperation: BinaryOperation {}

impl<Set: AssociativeBinaryOperation> AbstractAssociativeBinaryOperation<Set> for Set {}

// Commutativity
// ∀ a,b ∈ S, a∙b = b∙a

#[marker]
pub trait AbstractCommutativeBinaryOperation<Set: BaseSet>: AbstractBinaryOperation<Set> {}

#[marker]
pub trait CommutativeBinaryOperation: BinaryOperation {}

impl<Set: CommutativeBinaryOperation> AbstractCommutativeBinaryOperation<Set> for Set {}

// TODO: Macro that implements binary operation trait for types that satisfy std::Mul<Output = Self> trait

//==================//
// Identity Element //
//==================//

// Identity element for a binary operation
// ∃! e ∈ S, ∀ a ∈ S, a∙e = e∙a = a

pub trait AbstractIdentityElement<Set: BaseSet>: AbstractBinaryOperation<Set> {
    fn id() -> Set;
}

pub trait IdentityElement: BinaryOperation {
    fn id() -> Self;
    fn set_id(&mut self) {
        *self = Self::id();
    }
    fn is_id(&self) -> bool {
        self == &Self::id()
    }
}

impl<Set: IdentityElement> AbstractIdentityElement<Set> for Set {
    fn id() -> Set {
        Set::id()
    }
}

//=================//
// Inverse Element //
//=================//

// Inverse element for a binary operation
// ∀ a ∈ S, ∃! a' ∈ S, a∙a' = a'∙a = e

pub trait AbstractInverseElement<Set: BaseSet>:
    AbstractBinaryOperation<Set> + AbstractIdentityElement<Set>
{
    fn inv(a: Set) -> Set;
}

pub trait InverseElement: BinaryOperation + IdentityElement {
    fn inv(&self) -> Self;
}

impl<Set: InverseElement> AbstractInverseElement<Set> for Set {
    fn inv(a: Set) -> Set {
        a.inv()
    }
}

//============//
// Structures //
//============//

//=======//
// Magma //
//=======//

// Magma: set closed under a binary operation

pub trait AbstractMagma<Set: BaseSet> = AbstractBinaryOperation<Set>;

pub trait Magma = BinaryOperation;

//===========//
// Semigroup //
//===========//

// Semigroup: associative magma

pub trait AbstractSemigroup<Set: BaseSet> =
    AbstractMagma<Set> + AbstractAssociativeBinaryOperation<Set>;

pub trait Semigroup = Magma + AssociativeBinaryOperation;

//========//
// Monoid //
//========//

// Monoid: semigroup with identity element

pub trait AbstractMonoid<Set: BaseSet> = AbstractSemigroup<Set> + AbstractIdentityElement<Set>;

pub trait Monoid = Semigroup + IdentityElement;

//====================//
// Commutative monoid //
//====================//

// Commutative monoid: monoid whose binary operation is also commutative

pub trait AbstractCommutativeMonoid<Set: BaseSet> =
    AbstractMonoid<Set> + AbstractCommutativeBinaryOperation<Set>;

pub trait CommutativeMonoid = Monoid + CommutativeBinaryOperation;

//=======//
// Group //
//=======//

// Group: monoid for which every element has an inverse

pub trait AbstractGroup<Set: BaseSet> = AbstractMonoid<Set> + AbstractInverseElement<Set>;

pub trait Group = Monoid + InverseElement;

//===============//
// Abelian Group //
//===============//

// Abelian group: group whose binary operation is commutative

pub trait AbstractAbelianGroup<Set: BaseSet> =
    AbstractGroup<Set> + AbstractCommutativeBinaryOperation<Set>;

pub trait AbelianGroup = Group + CommutativeBinaryOperation;

pub mod instances {
    //! Example implementations of groupoids
    //!
    // TODO: Implement cool examples of groupoids.

    use super::*;

    // Cyclic group of order 2 (Z_2)
    // U∙A = A∙U = A
    // A ∙A = U
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub enum Z2 {
        U,
        A,
    }

    impl BinaryOperation for Z2 {
        fn op(&self, other: &Self) -> Self {
            match self {
                Z2::U => match other {
                    Z2::U => Z2::U,
                    Z2::A => Z2::A,
                },
                Z2::A => match other {
                    Z2::U => Z2::A,
                    Z2::A => Z2::U,
                },
            }
        }
    }

    impl AssociativeBinaryOperation for Z2 {}
    impl CommutativeBinaryOperation for Z2 {}

    impl IdentityElement for Z2 {
        fn id() -> Self {
            Z2::U
        }
    }

    impl InverseElement for Z2 {
        fn inv(&self) -> Self {
            match self {
                Z2::U => Z2::U,
                Z2::A => Z2::A,
            }
        }
    }
}

// Tests
#[cfg(test)]
mod tests {

    // TODO: Functions that test the instances defined above. Is there a way to auto test all of the properties?

    use super::{instances::*, AbelianGroup};

    // Test types
    fn is_abelian_group<Set: AbelianGroup>() {}

    #[test]
    fn groups() {
        //let u = Z2::U;
        //let a = Z2::A;

        is_abelian_group::<Z2>();
        // Test binary operation
        // Test unity
        // Test inverses
    }
}
