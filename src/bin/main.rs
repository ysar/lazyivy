use lazyivy::RungeKutta;
use ndarray::{Array, Array1};
 
 
fn brusselator(_t: &f64, y: &Array1<f64>) -> Array1<f64> {
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