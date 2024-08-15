//! # lazyivy
//!
//! lazyivy is a Rust crate that provides tools to solve initial value problems of the form
//! `dY/dt = F(t, Y)` using Runge-Kutta methods.
//!
//! The algorithms are implemented using the struct [`RungeKutta`] that implements
//! [`Iterator`]. The following Runge-Kutta methods are implemented currently, and
//! more will be added in the near future.  
//! - Euler 1
//! - Ralston 2
//! - Huen-Euler 2(1)
//! - Bogacki-Shampine 3(2)
//! - Fehlberg 4(5)
//! - Dormand-Prince 5(4)
//!
//! Where `p` is the order of the method and `p*` is the order of the error estimator step.
//!
//! ## Lazy integration  
//!
//! [`RungeKutta`] implements the `Iterator` trait. Each `.next()` call
//! advances the iteration to the next Runge-Kutta *step* and returns a tuple `(t, y)`, where `t`
//! is the dependent variable and `y` is `Array1<f64>`, which can be used to solve systems of ODEs.
//!
//! Note that each Runge-Kutta *step* contains `s` number of internal *stages*. Using lazyivy,
//! there is no way at present to access the integration values for these inner stages. The `next()`
//!  call returns the final result for each step, summed over all stages.
//!
//! The lazy implementation of Runge-Kutta means that you can consume the iterator in different
//! ways. For e.g., you can use `.last()` to keep only the final result, `.collect()` to gather the
//! state at all steps, `.map()` to chain the iterator with another, etc. You may also choose to
//! use it in a `for` loop and implement you own logic for modifying the step-size or customizing
//! the stop condition.
//!
//! ## Usage:
//!
//! After adding lazyivy to `Cargo.toml`, create an initial value problem using
//! the various `new_*` methods. Here is an example
//! showing how to solve the [Brusselator](https://en.wikipedia.org/wiki/Brusselator).  
//! ```rust
//! use lazyivy::RungeKutta;
//! use ndarray::{Array, ArrayView1, ArrayViewMut1};
//!  
//!  
//! fn brusselator(t: &f64, y: ArrayView1<f64>, mut result: ArrayViewMut1<f64>) {
//!     result[0] = 1. + y[0].powi(2) * y[1] - 4. * y[0];
//!     result[1] = 3. * y[0] - y[0].powi(2) * y[1];
//! }
//!  
//! fn main() -> Result<(), String> {
//!     let t0: f64 = 0.;
//!     let y0 = Array::from_vec(vec![1.5, 3.]);
//!     let absolute_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);
//!     let relative_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);
//!  
//!     // Instantiate a integrator for an ODE system with adaptive step-size Runge-Kutta.
//!  
//!     let mut integrator = RungeKutta::builder(brusselator, |t, _| *t > 40.)
//!         .initial_condition(t0, y0)
//!         .initial_step_size(0.025)
//!         .method("fehlberg", true)
//!         .tolerances(absolute_tol, relative_tol)
//!         .set_max_step_size(0.25)
//!         .build()?;
//!  
//!     // For adaptive algorithms, you can use this to improve the initial guess for the step size.
//!     integrator.set_step_size(&integrator.guess_initial_step());
//!  
//!     // Perform the iterations and print each state.
//!     for item in integrator {
//!         println!("{:?}", item)   // Prints the tuple (t, array[y1, y2]) at each iteration
//!     }
//!
//!     Ok(())
//! }
//! ```

#![warn(missing_docs)]

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
/// By convention the math community uses 1-based indexing to mark the a_ij matrix. Note that we
/// have used zero-based indexing above and in the implementation.  
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
