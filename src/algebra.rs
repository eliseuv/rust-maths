//! Simple representation of ususal algebraic structures using Rust's trait system
//!
//! `AbstractStrucuture<Set>` traits defines properties of a given algebraic structure over a given Set.
//! `Structure` traits represent trait objects equiped with a given algebraic structure. Here we call the set itself the algebraic structure.
//!
//! Every trait object also implements its corresponding abstract trait.
//! When an algebraic structure is built over another by imposing a property such associativity or commutativity, it is represented as a marker trait.
//!
//! Algebraic trait declaration goes as follows:
//! - Abstract trait declaration
//! - Trait object declaration
//! - Trait object implementation of abstract trait
//! - Implementation of trait object for objects that already implement certain traits
//! - Implementation of the previous trait, if it is a marker trait (This avoids having to manually implement marker traits on algebraic structures instances)

pub use std::ops::{Add, AddAssign, Mul, MulAssign, Neg};

//==========//
// Base Set //
//==========//

#[marker]
pub trait BaseSet: Sized + PartialEq {}
impl<Set> BaseSet for Set where Set: Sized + PartialEq {}

// Binary operation defined over a given set
// op: S ⨉ S → S
//     op(a,b) ↦ a·b
// ∀ a,b ∈ S, a∙b ∈ S
pub trait BinaryOperation<Set> {
    fn op(lhs: Set, rhs: Set) -> Set;
}

pub mod groupoid {
    //! Gropoids are group-like algebraic structures.
    //! They are algebraic structures built from a set closed under a binary operation.
    //!
    //! Groupoid algebraic structures:
    //! - Magma: set closed under a binary operation
    //! - Semigroup: associative magma
    //! - Monoid: semigroup with unit element
    //! - Group: monoid with an inverse for every element
    //! - Abelian group: commutative group

    //=======//
    // Magma //
    //=======//

    // Magma (or Groupoid): set closed under a binary operation
    // ∃ op: S ⨉ S → S
    //    ·: (a,b) ↦ a·b
    // ∀ a,b ∈ S, a∙b ∈ S
    pub trait AbstractMagma<Set> {
        fn op(lhs: Set, rhs: Set) -> Set;
    }

    // Magma trait object
    pub trait Magma: Sized + super::Mul<Output = Self> + super::MulAssign {}

    // Magma trait object also implements magma abstract trait
    impl<Set: Magma> AbstractMagma<Set> for Set {
        fn op(lhs: Self, rhs: Self) -> Self {
            lhs * rhs
        }
    }

    // The magma trait is automatically satisfied for every sized object that has multiplication operator defined
    impl<Set> Magma for Set where Set: Sized + super::Mul<Output = Self> + super::MulAssign {}

    //===========//
    // Semigroup //
    //===========//

    // Semigroup: associative magma
    // ∀ a,b,c ∈ S, (a∙b)∙c = a∙(b∙c)
    #[marker]
    pub trait AbstractSemigroup<Set>: AbstractMagma<Set> {}

    // Semigroup trait object
    #[marker]
    pub trait Semigroup: Magma {}

    // Semigroup trait object also implements abstract trait
    impl<Set: Semigroup> AbstractSemigroup<Set> for Set {}

    //========//
    // Monoid //
    //========//

    // Monoid: semigroup with identity element
    // ∃! e ∈ S, ∀ a ∈ S, a∙e = e∙a = a
    pub trait AbstractMonoid<Set>: AbstractSemigroup<Set> {
        // Get identity
        fn id() -> Set;
    }

    // Monoid trait object
    pub trait Monoid: Semigroup + AbstractMonoid<Self> {
        // Get identity
        fn get_id() -> Self;
        // Set to identity
        fn set_id(&mut self) {
            *self = Self::id();
        }
        // Test if is identity
        fn is_id(&self) -> bool;
    }

    // Monoid trait object also implements monoid abstract trait
    impl<Set: Monoid> AbstractMonoid<Set> for Set {
        fn id() -> Set {
            Monoid::get_id()
        }
    }

    // Implement semigroup for monoids
    impl<T, Set> AbstractSemigroup<Set> for T where T: AbstractMonoid<Set> {}
    impl<Set: Monoid> Semigroup for Set {}

    //====================//
    // Commutative monoid //
    //====================//

    //=======//
    // Group //
    //=======//

    // Group is a monoid for which every element has an inverse
    // ∀ a ∈ S, ∃! a' ∈ S, a∙a' = a'∙a = e
    pub trait AbstractGroup<Set> {
        fn inv(a: &Set) -> Set;
    }

    // Group trait object
    pub trait Group: Monoid {
        fn inv(&self) -> Self;
    }

    // Group trait object also implements abstract trait
    impl<Set: Group> AbstractGroup<Set> for Set {
        fn inv(a: &Set) -> Set {
            a.inv()
        }
    }

    //===============//
    // Abelian Group //
    //===============//

    // Abelian group is a group whose binary operation is commutative
    // ∀ a,b ∈ S, a∙b = b∙a
    #[marker]
    pub trait AbstractAbelianGroup<Set>: AbstractGroup<Set> {}

    // Abelian group trait object
    #[marker]
    pub trait AbelianGroup: Group {}

    // Abelian group trait object also implements abstract trait
    impl<Set: AbelianGroup> AbstractAbelianGroup<Set> for Set {}

    pub mod instances {
        //! Example implementations of groupoids
        //!

        use super::{super::Mul, super::MulAssign, AbelianGroup, Group, Monoid};

        // Cyclic group of order 2 (Z_2)
        // U ∙ A = A ∙ U = A
        // A ∙ A = U
        #[derive(Debug, Clone, Copy, PartialEq)]
        pub enum Z2 {
            U,
            A,
        }

        impl Mul<Self> for Z2 {
            type Output = Self;

            fn mul(self, rhs: Self) -> Self::Output {
                match self {
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

        impl MulAssign for Z2 {
            fn mul_assign(&mut self, rhs: Self) {
                *self = self.mul(rhs)
            }
        }

        impl Monoid for Z2 {
            fn get_id() -> Self {
                Z2::U
            }
            fn is_id(&self) -> bool {
                self == &Z2::U
            }
        }

        impl Group for Z2 {
            fn inv(&self) -> Self {
                match self {
                    Z2::U => Z2::U,
                    Z2::A => Z2::A,
                }
            }
        }

        impl AbelianGroup for Z2 {}
    }

    // Tests
    #[cfg(test)]
    mod tests {

        use super::instances::*;
        use super::*;

        #[test]
        fn groups() {
            let u = Z2::U;
            let a = Z2::A;
            // Test binary operation
            assert_eq!(u * u, u);
            assert_eq!(u * a, a);
            assert_eq!(a * u, a);
            assert_eq!(a * a, u);
            // Test unity
            assert_eq!(Z2::id(), u);
            // Test inverses
            assert_eq!(u * u.inv(), Z2::id());
            assert_eq!(a * a.inv(), Z2::id());
        }
    }
}

pub mod ringoid {
    //! Ringoids are ring-like algebraic structures.
    //! They are algebraic structures built from a set closed under two binary operations, often called addition and multiplication, with multiplication being distributive over addition.
    //!
    //! Ringoid algebraic structures:
    //! - Ringoid: set with two binary operations, often called addition and multiplication, with multiplication being distributive over addition
    //! - Semiring: ringoid such that is also a monoid under each binary operation (the multiplicative identity is usually called 1 and the additive identity is usually called 0)
    //! - Near ring: semiring whose addition monoid is a group
    //! - Ring: semiring whose addition monoid is an abelain group
    //! - Commutative ring: ring whose multiplication operation is commutative
    //! - Field: Commutative ring which contains a multiplicative inverse for every nonzero element (i.e., except for the identity of the addition monoid)

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
    pub trait Ringoid:
        Sized
        + super::Mul<Output = Self>
        + super::MulAssign
        + super::Add<Output = Self>
        + super::AddAssign
    {
    }

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
        Set: Sized
            + super::Mul<Output = Self>
            + super::MulAssign
            + super::Add<Output = Self>
            + super::AddAssign
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
    pub trait NearRing: Ringoid + AbstractSemiring<Self> + super::Neg<Output = Self> {
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
}
