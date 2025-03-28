use criterion::{criterion_group, criterion_main, Criterion};
use lazyivy::{RungeKutta, RungeKuttaMethod};
use ndarray::{array, ArrayView1, ArrayViewMut1};
use paste::paste;

fn brusselator(_t: &f64, y: ArrayView1<f64>, mut result: ArrayViewMut1<f64>) {
    result[0] = 1. + y[0].powi(2) * y[1] - 4. * y[0];
    result[1] = 3. * y[0] - y[0].powi(2) * y[1];
}

macro_rules! bench_rungekutta_func {
    ($($method:ident), *) => {
        $(
            paste! {

                fn [< bench_ $method:lower >]() {
                    let t0: f64 = 0.;
                    let y0 = array![1.5, 3.];
                    let absolute_tol = array![1.0e-4, 1.0e-4];
                    let relative_tol = array![1.0e-4, 1.0e-4];

                    let integrator = RungeKutta::builder(brusselator, |t, _| *t > 200.)
                        .initial_condition(t0, y0)
                        .initial_step_size(0.025)
                        .method(RungeKuttaMethod::$method, false)
                        .tolerances(absolute_tol, relative_tol)
                        // .set_max_step_size(0.25)
                        .build()
                        .unwrap();

                    // integrator.set_step_size(&integrator.guess_initial_step());

                    // Perform the iterations
                    for _item in integrator { }
                }
            }
        )*
    };
}

bench_rungekutta_func!(Euler, Ralston, HuenEuler, BogackiShampine, Fehlberg, DormandPrince);

fn criterion_benchmark_rungekutta(c: &mut Criterion) {
    let mut group = c.benchmark_group("Runge-Kutta");
    group.bench_function("Euler", |b| b.iter(|| bench_euler()));
    group.bench_function("Ralston", |b| b.iter(|| bench_ralston()));
    group.bench_function("HuenEuler", |b| b.iter(|| bench_hueneuler()));
    group.bench_function("Bogacki-Shampine", |b| b.iter(|| bench_bogackishampine()));
    group.bench_function("Fehlberg", |b| b.iter(|| bench_fehlberg()));
    group.bench_function("Dormand-Prince", |b| b.iter(|| bench_dormandprince()));
    group.finish();
}

criterion_group!(benches, criterion_benchmark_rungekutta);
criterion_main!(benches);
