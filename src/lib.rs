//! # lazyivy
//!
//! lazyivy is a crate that provides tools to solve initial value problems of
//! the form `dY/dt = F(t, y)` using Runge-Kutta methods.  
//! An initial value problem can be defined by instantiating a struct of type -
//! - [`RungeKutta`](explicit_single::RungeKutta),
//!   for using Runge-Kutta methods with a fixed step size.
//!- [`RungeKuttaAdaptive`](explicit_single::RungeKuttaAdaptive),
//!   for using Runge-Kutta methods with an adaptive step size.

#![warn(missing_docs)]

/// Auxiliary methods
#[doc(hidden)]
pub mod aux;

/// The Butcher table contains the `a_ij` matrix and the `c_i` and `b_j` vector
/// coefficients that specify a particular Runge-Kutta method.  
///
/// It is usually written like this -
/// ```
/// ---------------------------------------------
/// c0     | -
/// c1     | a10, -
/// c2     | a20, a21, -
/// c3     | a30, a31, a32, -
/// ...
/// c{s-1} | a{s-1}0, a{s-1}1, ... a{s-1}{s-2}, -
/// ---------------------------------------------
///        | b0, b1, b2, ..., b{s-1}
/// ---------------------------------------------
/// ```
/// By convention the math community uses 1-based indexing to mark the a_ij
/// matrix. Note that we have used zero-based indexing above and in the
/// implementation.  
///
/// For adaptive Runge-Kutta methods the lower-order coefficients are also
/// provided for estimating the error
/// ```
///        | b20, b21, b22, ..., b2{s-1}
/// ---------------------------------------------
/// ```
/// All coefficients are stored as 1D arrays. For explicit Runge-Kutta methods,
/// `a_ij` is a lower triangular matrix and is stored as a one-dimension `array`
/// with indexing logic `idx = i(i-1)/2 + j`.
///
/// In the Runge-Kutta iteration, we initialize `k(i=0)` with the state at the
/// start of the iteration. We start the stages from i=1, since `a00` is
/// undefined for any explicit Runge-Kutta Method.
pub mod tables;

/// Runge-Kutta methods for a single variable ODE.
pub mod explicit_single;

/// Runge-Kutta methods for a system of ODEs, i.e. with multiple variables.
pub mod explicit_system;
