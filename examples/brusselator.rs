use lazyivy::RungeKutta;
use ndarray::{array, ArrayView1, ArrayViewMut1};
use std::fs::File;
use std::io::Write;

fn brusselator(_t: &f64, y: ArrayView1<f64>, mut result: ArrayViewMut1<f64>) {
    result[0] = 1. + y[0].powi(2) * y[1] - 4. * y[0];
    result[1] = 3. * y[0] - y[0].powi(2) * y[1];
}

fn main() {
    let t0: f64 = 0.;
    let y0 = array![1.5, 3.];
    let absolute_tol = array![1.0e-4, 1.0e-4];
    let relative_tol = array![1.0e-4, 1.0e-4];

    let integrator = RungeKutta::builder(brusselator, |t, _| *t > 40.)
        .initial_condition(t0, y0)
        .initial_step_size(0.025)
        .method("dormandprince", true)
        .tolerances(absolute_tol, relative_tol)
        .set_max_step_size(0.25)
        .build()
        .unwrap();

    let mut buffer = File::create("data/brusselator.csv").unwrap();

    for item in integrator {
        // writeln!(
        //     &mut buffer,
        //     "{:.4}, {:.4}, {:.4}",
        //     item.0, item.1[0], item.1[1]
        // )
        // .unwrap();
    }

}
