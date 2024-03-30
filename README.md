# lazyivy

[![Crate](https://img.shields.io/crates/v/lazyivy)](https://crates.io/crates/lazyivy)
[![Build](https://github.com/ysar/lazyivy/actions/workflows/build.yml/badge.svg)](https://github.com/ysar/lazyivy/actions/workflows/build.yml)
[![Documentation](https://img.shields.io/docsrs/lazyivy/latest)](https://docs.rs/lazyivy/latest/lazyivy/)

lazyivy is a Rust crate that provides tools to solve initial value problems of 
the form `dY/dt = F(t, y)` using Runge-Kutta methods. 

The algorithms are implemented using the struct `RungeKutta`, that implements 
`Iterator`. The following Runge-Kutta methods are implemented currently, and 
more will be added in the near future.  
- Euler 1
- Ralston 2
- Huen-Euler 2(1)
- Bogacki-Shampine 3(2)
- Fehlberg 4(5)
- Dormand-Prince 5(4)

Where `p` is the order of the method and `(p*)` is the order of the embedded 
error estimator, if it is present.

## Lazy integration
`RungeKutta` implements the `Iterator` trait. Each `.next()` call advances the 
iteration to the next Runge-Kutta *step* and returns a tuple `(t, y)`, where 
`t` is the dependent variable and `y` is `Array1<f64>`. 

Note that each Runge-Kutta *step* contains `s` number of internal *stages*. 
Using lazyivy, there is no way at present to access the integration values for 
these inner stages. The `.next()` call returns the final result for each step, 
summed over all stages.

The lazy implementation of Runge-Kutta means that you can consume the iterator 
in different ways. For e.g., you can use `.last()` to keep only the final 
result, `.collect()` to gather the state at all steps, `.map()` to chain the 
iterator with another, etc. You may also choose to use it in a `for` loop and 
implement you own logic for modifying the step-size or customizing the stop 
condition.

**API is unstable. It is active and under development.**

## Usage: 

After adding lazyivy to `Cargo.toml`, create an initial value problem using 
the various `new_*` methods. Here is an example 
showing how to solve the [Brusselator](https://en.wikipedia.org/wiki/Brusselator). 

```math 
\frac{d}{dt} \left[ \begin{array}{c}
 y_1 \\ y_2 \end{array}\right] = \left[\begin{array}{c}1 - y_1^2 y_2 - 4 y_1 
 \\ 3y_1 - y_1^2 y_2 \end{array}\right]
```
```rust
use lazyivy::RungeKutta;
use ndarray::{Array, Array1};
 
 
fn brusselator(t: &f64, y: &Array1<f64>) -> Array1<f64> {
    Array::from_vec(vec![
        1. + y[0].powi(2) * y[1] - 4. * y[0],
        3. * y[0] - y[0].powi(2) * y[1],
    ])
}
 
fn main() {
    let t0: f64 = 0.;
    let y0 = Array::from_vec(vec![1.5, 3.]);
    let absolute_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);
    let relative_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);
 
    // Instantiate a integrator for an ODE system with adaptive step-size 
    // Runge-Kutta.
 
    let mut integrator = RungeKutta::new_fehlberg(
        t0,              // Initial condition - time t0
        y0,              // Initial condition - Initial condition [y1, y2] @ t0
        brusselator,     // Evaluation function
        |t, _| t > &20., // Predicate that determines stop condition
        0.025,           // Initial step size
        relative_tol,    // Relative tolerance for error estimation
        absolute_tol,    // Absolute tolerance for error estimation
        true,            // Use adaptive step-size
    );
 
    // For adaptive algorithms, you can use this to improve the initial guess 
    // for the step size.
    integrator.set_step_size(&integrator.guess_initial_step());
 
    // Perform the iterations and print each state.
    for item in integrator {
        println!("{:?}", item)   // Prints (t, array[y1, y2]) for each step.
    }
}
```
The result when plotted looks like this - 
![Brusselator](https://raw.githubusercontent.com/ysar/lazyivy/main/examples/brusselator_adaptive.png)

## To-do list:
- [ ] Add more Runge-Kutta methods.
- [ ] Improve tests.
- [ ] Benchmark.
- [x] Move allocations out of `next` and into a separate struct.
- [ ] Add more examples, e.g. Lorentz attractor.
- [ ] 