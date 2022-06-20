//! Simple topological structures

use crate::sets::FiniteSet;
use num_traits::real::Real;

/// Finite metric space
///
/// Metric
///
/// The metric function is defined by any pair of set elements and returns a positive real:
///
///     d: S ⨉ S → ℝ⁺
///
/// such that ∀ x,y,z ∈ S:
///
/// - Identity of indiscernables
///     d(x,y) = 0  ⇔  x = y
///
/// - Symmetry
///     d(x,y) = d(y,x)
///
/// - Triangle inequality
///     d(x,z) ≤ d(x,y) + d(y,z)

pub trait AbstractFiniteMetricSpace<Set, PositiveReal>
where
    Set: FiniteSet,
    PositiveReal: Real,
{
    fn metric(lhs: &Set, rhs: &Set) -> PositiveReal;
}

pub trait FiniteMetricSpace<PositiveReal: Real>: FiniteSet {
    fn metric(&self, other: &Self) -> PositiveReal;
}

impl<Set, PositiveReal> AbstractFiniteMetricSpace<Set, PositiveReal> for Set
where
    Set: FiniteMetricSpace<PositiveReal>,
    PositiveReal: Real,
{
    fn metric(lhs: &Set, rhs: &Set) -> PositiveReal {
        lhs.metric(rhs)
    }
}
