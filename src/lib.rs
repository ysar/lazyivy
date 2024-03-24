//! # lazyivy
//!
//! lazyivy is a crate that provides tools to solve initial value problems of the form
//! `dY/dt = F(t, y)` using Runge-Kutta methods.
//!
//! Fixed step-size Runge-Kutta algorithms are implemented using the the structs,
//! - [`RungeKutta`](explicit_single::RungeKutta) to solve an ODE of one variable.
//! - [`RungeKuttaSystem`](explicit_system::RungeKuttaSystem) to solve a system of ODEs.  
//!
//! These include the following RK methods - `Euler`, `Ralston`
//!
//! Also implemented are embedded Runge-Kutta methods that consider an adaptive step-size based on
//! the local error estimate from a lower order method. These are implemented using the structs,
//! - [`RungeKuttaAdaptive`](explicit_single::RungeKuttaAdaptive) to solve an ODE of one variable
//!   using an adaptive step size.
//! - [`RungeKuttaSystemAdaptive`](explicit_system::RungeKuttaSystemAdaptive) to solve a system of
//!   ODEs with an adaptive step size.
//!
//! These include the following RK methods - `Fehlberg`, `Hueneuler`, `DormandPrince`.
//!
//! # Usage:
//! Instantiate a struct object and call the `new_` methods corresponding to the Runge-Kutta method
//! of your choice. The `new_` method calls for `RungeKutta` and `RungeKuttaSystem` take the
//! following form.
//! ```rust
//! let euler_integrator = RungeKutta::new_euler(t0, y0, f, p, step);
//! let euler_integrator_system = RungeKuttaSystem::new_euler(t0, y0, f, p, step);
//! ```
//! Where you can replace `_euler` with any other Runge-Kutta method that is implemented for fixed
//! step-size.

#![warn(missing_docs)]

/// Auxiliary methods
#[doc(hidden)]
pub mod aux;

/// The Butcher table contains the `a_ij` matrix and the `c_i` and `b_j` vector coefficients that
/// specify a particular Runge-Kutta method.  
///
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
/// By convention the math community uses 1-based indexing to mark the a_ij matrix. Note that we
/// have used zero-based indexing above and in the implementation.  
///
/// For adaptive Runge-Kutta methods the lower-order coefficients are also provided for estimating
/// the error.
/// ```
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

/// Runge-Kutta methods for a single variable ODE.
pub mod explicit_single;

/// Runge-Kutta methods for a system of ODEs, i.e. with multiple variables.
pub mod explicit_system;
