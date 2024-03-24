# lazyivy

lazyivy is a Rust crate that provides tools to solve initial value problems of the form
`dY/dt = F(t, y)` using Runge-Kutta methods. 

Fixed step-size Runge-Kutta algorithms are implemented using the the structs,
- [`RungeKutta`](explicit_single::RungeKutta) to solve an ODE of one variable, and
- [`RungeKuttaSystem`](explicit_system::RungeKuttaSystem) to solve a system of ODEs.  

These include the following RK methods - Euler, Ralston

Also implemented are embedded Runge-Kutta methods that consider an adaptive step-size based on
the local error estimate from a lower order method. These are implemented using the structs,
- [`RungeKuttaAdaptive`](explicit_single::RungeKuttaAdaptive) to solve an ODE of one variable 
  using an adaptive step size.
- [`RungeKuttaSystemAdaptive`](explicit_system::RungeKuttaSystemAdaptive) to solve a system of ODEs 
  with an adaptive step size.

These include the following RK methods of order `p(p*)` - Fehlberg 4(5), Hueneuler, Dormand-Prince 4(5).  

Where `p` is the order of the method and `p*` is the order of the error estimator step.

All integration structs implement the `Iterator` trait. Each `.next()` call advances the iteration
to the next Runge-Kutta *step* and returns a tuple `(t, y)`, where `t` is the dependent variable and 
`y` can be either `f64`, as in the case of `RungeKutta` or `Array1<f64>`, in the case of 
`RungeKuttaSystem`. 

Note that each Runge-Kutta *step* contains `s` number of internal *stages*. Using lazyivy, there is 
no way at present to access the integration values for these inner stages. The `next()` call returns 
the final result for each step, summed over all stages.

The lazy implementation of Runge-Kutta means that you have a range of options to consume the 
iterator. For e.g., you can use `.last()` to keep only the final result, `.collect()` to gather the
state at all steps, `.map()` to chain the iterator with another, etc. You may also choose to use it 
in a `for` loop and implement you own logic for modifying the step-size or customizing the stop 
condition.

## Usage: 

After adding lazyivy to `Cargo.toml`, create an initial value problem using 
the various `new_*` methods. Here is an example 
showing how to solve the [Brusselator](https://en.wikipedia.org/wiki/Brusselator). 
```math
\frac{d}{dt} \left[ \begin{array}{c}
 y_1 \\ y_2 \end{array}\right] = \left[\begin{array}{c}1 - y_1^2 y_2 - 4 y_1 \\ 3y_1 - y_1^2 y_2 \end{array}\right]
```
```rust
use lazyivy::RungeKuttaSystemAdaptive;
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
 
    // Instantiate a integrator for an ODE system with adaptive step-size Runge-Kutta.
 
    let mut integrator = RungeKuttaSystemAdaptive::new_fehlberg(
        t0,                   // Initial condition - time
        y0,                   // Initial condition - Brusselator variables in Array[y1, y2]
        brusselator,          // Evaluation function
        |t, _| t > &20.,      // Predicate that determines stop condition
        0.025,                // Initial step size
        relative_tol,         // Relative tolerance for error estimation
        absolute_tol,         // Absolute tolerance for error estimation
    );
 
    // For adaptive algorithms, you can use this to improve the initial guess for the step size.
    integrator.h = integrator.guess_initial_step();
 
    // Perform the iterations and print each state.
    for item in integrator {
        println!("{:?}", item)   // Prints the tuple (t, array[y1, y2]) at each iteration 
    }
}
```
The result when plotted looks like this - 
![Brusselator](examples/brusselator_adaptive.png)