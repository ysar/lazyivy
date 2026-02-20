use lazyivy::{RungeKutta, RungeKuttaMethod};
use ndarray::{array, Array1, ArrayView1, ArrayViewMut1};
use std::fs::File;
use std::io::Write;


fn lorentz_attractor(
    x: f64,
    y: f64,
    z: f64,
    sigma: f64,
    beta: f64,
    rho: f64
    ) -> Array1<f64> {
    array![
        sigma * (y - x),
        x * (rho - z) - y,
        x * y - beta * z,
    ]
}

fn main() {
    let t0: f64 = 0.;
    let y0 = array![2., 1., 1.];
    let absolute_tol = array![1.0e-4, 1.0e-4, 1.0e-4];
    let relative_tol = array![1.0e-4, 1.0e-4, 1.0e-4];

    let sigma = 10.;
    let beta = 8. / 3.;
    let rho = 28.;

    let eval_closure = |_t: &f64, y: ArrayView1<f64>, mut result: ArrayViewMut1<f64>| {
        // Closure captures the environment and wraps the function signature
        result.assign(&(lorentz_attractor(y[0], y[1], y[2], sigma, beta, rho)));
    };

    let integrator = RungeKutta::builder(eval_closure, |t, _| *t > 30.)
        .initial_condition(t0, y0)
        .initial_step_size(0.001)
        .method(RungeKuttaMethod::DormandPrince, true)
        .tolerances(absolute_tol, relative_tol)
        .set_max_step_size(0.01)
        .build()
        .unwrap();

    let mut buffer = File::create("data/lorenz_attractor.csv").unwrap();

    for item in integrator {
        writeln!(
            &mut buffer,
            "{:.4}, {:.4}, {:.4}, {:.4}",
            item.0, item.1[0], item.1[1], item.1[2],
        )
        .unwrap();
    }
}
