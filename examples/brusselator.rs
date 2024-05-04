use lazyivy::RungeKutta;
use ndarray::{array, Array1};
use std::fs::File;
use std::io::Write;

fn brusselator(_t: &f64, y: &Array1<f64>) -> Array1<f64> {
    array![
        1. + y[0].powi(2) * y[1] - 4. * y[0],
        3. * y[0] - y[0].powi(2) * y[1],
    ]
}

fn main() -> Result<(), String> {
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
        .build()?;

    let mut buffer = File::create("data/brusselator.csv").unwrap();

    for item in integrator {
        writeln!(
            &mut buffer,
            "{:.4}, {:.4}, {:.4}",
            item.0, item.1[0], item.1[1]
        )
        .unwrap();
    }

    Ok(())
}
