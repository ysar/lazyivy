# lazyivy

lazyivy is a Rust library for solving initial value problems (IVPs) 
by using Runge-Kutta integration implemented using iterators. It provides a 
struct called `RungeKuttaIntegrator` that implements the `Iterator` trait. 

## Usage: 

After adding lazyivy to `Cargo.toml`, create an initial value problem using 
the various `new_*` methods for `RungeKuttaIntegrator`. Here is an example 
showing how to integrate `dy/dx = y` using the Euler method. 

```rust

let t0 = 0.;
let y0 = 0.;
let h = 1.;
let integrator = RungeKuttaIntegrator::new_euler(t0, y0, |t, y| y, h);
let (t_result, y_result) = integrator.last(); 
    // Consumes iterator and returns last value 
```