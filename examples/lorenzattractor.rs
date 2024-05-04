use lazyivy::RungeKutta;
use ndarray::{array, Array1};
use std::fs::File;
use std::io::Write;

fn lorentz_attractor(_t: &f64, y: &Array1<f64>) -> Array1<f64> {
    array![
        10. * (y[1] - y[0]),
        y[0] * (28. - y[2]) - y[1],
        y[0] * y[1] - 8. / 3. * y[2],
    ]
}

fn main() -> Result<(), String> {
    let t0: f64 = 0.;
    let y0 = array![2., 1., 1.];
    let absolute_tol = array![1.0e-4, 1.0e-4, 1.0e-4];
    let relative_tol = array![1.0e-4, 1.0e-4, 1.0e-4];

    let integrator = RungeKutta::builder(lorentz_attractor, |t, _| *t > 30.)
        .initial_condition(t0, y0)
        .initial_step_size(0.001)
        .method("dormandprince", true)
        .tolerances(absolute_tol, relative_tol)
        .set_max_step_size(0.01)
        .build()?;

    let mut buffer = File::create("data/lorentz_attractor.csv").unwrap();

    for item in integrator {
        writeln!(
            &mut buffer,
            "{:.4}, {:.4}, {:.4}, {:.4}",
            item.0, item.1[0], item.1[1], item.1[2],
        )
        .unwrap();
    }

    Ok(())
}
