//! # lazyivy
//!
//! lazyivy is a crate that provides tools to solve initial value problems of
//! the form `dY/dt = F(t, y)` using Runge-Kutta methods.  
//! An initial value problem can be defined by instantiating a struct of type -
//! - [`RungeKutta`](explicit_float::RungeKutta),
//!   for using Runge-Kutta methods with a fixed step size.
//!- [`RungeKuttaAdaptive`](explicit_float::RungeKuttaAdaptive),
//!   for using Runge-Kutta methods with an adaptive step size.

pub mod butcher;
pub mod explicit_float;
