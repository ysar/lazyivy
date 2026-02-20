#![doc=include_str!("../README.md")]
#![warn(missing_docs)]

/// Contains the `BuilderError` error struct that is returned if the Runge-Kutta builder fails.
pub mod error;

/// The Butcher table contains the `a_ij` matrix and the `c_i` and `b_j` vector coefficients that
/// specify a particular Runge-Kutta method.
///
/// ```text
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
/// By convention the math community uses 1-based indexing to mark the a_ij matrix. Note - this lib
/// uses zero-based indexing above and in the implementation.
///
/// For adaptive Runge-Kutta methods the lower-order coefficients are also provided for estimating
/// the error.
/// ```text
///        | b20, b21, b22, ..., b2{s-1}
/// ---------------------------------------------
/// ```
/// All coefficients are stored as 1D arrays. For explicit Runge-Kutta methods, `a_ij` is a lower
/// triangular matrix and is stored as a one-dimension `array` with indexing logic
/// `idx = i(i-1)/2 + j`.
///
/// In the Runge-Kutta iteration, we initialize `k(i=0)` with the state at the start of the
/// iteration. We start the stages from i=1, since `a00` is undefined for any explicit
/// Runge-Kutta Method.
pub mod tables;

/// Explicit Runge-Kutta Methods
pub mod explicit;

// Re-export some useful structs
pub use crate::explicit::RungeKutta;
pub use crate::tables::RungeKuttaMethod;
