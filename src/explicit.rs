use crate::misc;
use crate::tables;
use ndarray::Array1;
use paste::paste;

macro_rules! impl_new_rungekutta {
    ($name:ident) => {
        paste! {
            #[doc="Instantiate an ODE integrator (single variable) which uses the " $name:camel " method with adaptive step-size."]
            pub fn [<new_ $name:lower >] (
                t: f64,
                y: Array1<f64>,
                f: F,
                predicate: P,
                h: f64,
                relative_tol: Array1<f64>,
                absolute_tol: Array1<f64>,
                do_adaptive: bool,
            ) -> Self {

                let num_variables = y.len();
                let zeros = Array1::<f64>::zeros(num_variables);
                let table = tables::[<$name:upper _BT>];

                // Pre-allocate some arrays used in the iterations. Cloning here is fine, we are
                // only doing it once at the beginning.
                let tmp_variables = _TempArrays {
                    k: vec![zeros.clone(); table.s],
                    y: zeros.clone(),
                    y_next: zeros.clone(),
                    y_err: zeros.clone(),
                    sum: zeros.clone()
                };

                // Instantiate and return RungeKutta instance.
                RungeKutta {
                    t,
                    y: y.clone(),
                    f,
                    predicate,
                    h,
                    relative_tol,
                    absolute_tol,
                    do_adaptive,
                    table,
                    _num_variables: num_variables,
                    _t: tmp_variables
                }
            }
        }
    }
}

/// Struct to store preallocated arrays
struct _TempArrays {
    k: Vec<Array1<f64>>,
    y: Array1<f64>,
    y_next: Array1<f64>,
    y_err: Array1<f64>,
    sum: Array1<f64>,
}

/// Struct for adaptive step size Runge-Kutta integration. Stores the integration state at every
/// iteration. Implements [`Iterator`].
///
/// ## Usage:
/// ```rust

/// ```
pub struct RungeKutta<'a, F, P>
where
    F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
    P: Fn(&f64, &Array1<f64>) -> bool,
{
    t: f64,
    y: Array1<f64>,
    f: F,
    predicate: P,
    h: f64,
    relative_tol: Array1<f64>,
    absolute_tol: Array1<f64>,
    do_adaptive: bool,
    table: tables::ButcherTableau<'a>,
    _num_variables: usize,
    _t: _TempArrays,
}

///
impl<F, P> RungeKutta<'_, F, P>
where
    F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
    P: Fn(&f64, &Array1<f64>) -> bool,
{
    impl_new_rungekutta!(Euler);
    impl_new_rungekutta!(Ralston);
    impl_new_rungekutta!(HuenEuler);
    impl_new_rungekutta!(BogackiShampine);
    impl_new_rungekutta!(Fehlberg);
    impl_new_rungekutta!(DormandPrince);

    /// Set step-size.
    pub fn set_step_size(&mut self, h: &f64) {
        self.h = *h;
    }

    /// Calculates the norm || (y0 - y1) || as defined in Harrier, Nørsett, Wanner.
    fn calc_error_norm(&self, y0: &Array1<f64>, y1: &Array1<f64>) -> f64 {
        let mut tolerance = Array1::<f64>::zeros(self._num_variables);
        for i in 0..self._num_variables {
            tolerance[i] =
                self.absolute_tol[i] + y0[i].abs().max(y1[i].abs()) * self.relative_tol[i];
        }
        (1. / self._num_variables as f64 * ((y0 - y1) / tolerance).map(|a| a * a).sum()).sqrt()
    }

    /// Provides a reasonable guess for the step-size at the first iteration.
    /// Follows the algorithm written in Harrier, Nørsett, Wanner.
    pub fn guess_initial_step(&self) -> f64 {
        let f0 = (self.f)(&self.t, &self.y);
        let zero_array = Array1::<f64>::zeros(self._num_variables);
        let d0 = self.calc_error_norm(&self.y, &zero_array);
        let d1 = self.calc_error_norm(&f0, &zero_array);

        let h0: f64 = if d0 < 1.0e-5 || d1 < 1.0e-5 {
            1.0e-6
        } else {
            0.01 * (d0 / d1)
        };

        let y1 = &self.y + h0 * &f0;
        let f1 = (self.f)(&(self.t + h0), &y1);

        let d2 = self.calc_error_norm(&f1, &f0) / h0;

        let h1: f64 = if d1.max(d2) < 1e-15 {
            (h0 * 1.0e-3).max(1.0e-6)
        } else {
            (0.01 / d1.max(d2)).powf(1. / (self.table.p + 1) as f64)
        };

        h1.min(100. * h0)
    }
}

impl<F, P> Iterator for RungeKutta<'_, F, P>
where
    F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
    P: Fn(&f64, &Array1<f64>) -> bool,
{
    type Item = (f64, Array1<f64>);

    fn next(&mut self) -> Option<Self::Item> {
        // Check stop condition at every iteration
        if (self.predicate)(&self.t, &self.y) {
            return None;
        }

        // RK Logic

        // Store k[i] values in a vector, initialize with the first value
        let f_now = (self.f)(&self.t, &self.y);

        // Re-fill the temporary arrays with zeros or fill values
        self._t.k.fill(f_now);
        self._t.y.fill(0.);
        self._t.y_next.fill(0.);
        self._t.y_err.fill(0.);

        #[allow(unused_assignments)]
        let mut t: f64 = 0.;

        #[allow(unused_assignments)]
        let mut indx: usize = 0;

        #[allow(unused_assignments)]
        let mut t_next: f64 = 0.;

        let mut h = self.h;

        // Set the error norm to be larger than 1. before first while condition check.
        let mut error_norm = 2.;

        // Scale factors for changing the step-size. Prevents overzealous changes.
        let scale_factor: f64 = 0.9;
        let min_scale_factor: f64 = 0.1;
        let max_scale_factor: f64 = 2.;

        while error_norm > 1. {
            self._t.sum.fill(0.);

            for i in 1..self.table.s {
                // Index into the lower triangular matrix.
                indx = i * (i - 1) / 2;

                // Advance time in the inner stage.
                t = self.t + self.table.c[i] * self.h;

                // Calculate `y` in the inner stage using k[0], k[1], ... k[stage - 1].
                self._t.y = &self.y
                    + self.h
                        * misc::sum_product_array(&self.table.a[indx..indx + 1], &self._t.k[..i]);

                // Calculate `k` for the next stage.
                self._t.k[i] = (self.f)(&t, &self._t.y);

                // Add the stage result and move to next stage.
                self._t.sum = &self._t.sum + self.table.b[i] * &self._t.k[i];
            }

            // Advance time and output in the outer stage.
            t_next = self.t + h;
            self._t.y_next = &self.y + h * self.table.b[0] * &self._t.k[0] + h * &self._t.sum;

            // If no adaptive step-size requested, we can break out of the loop now.
            if !self.do_adaptive {
                break;
            }

            // Calculate the result using the error estimator coefficients.
            self._t.y_err = &self.y + h * misc::sum_product_array(self.table.b2, &self._t.k);

            error_norm = self.calc_error_norm(&self._t.y_next, &self._t.y_err);

            // Scale step-size based on the norm of the error between the result and the error
            // estimator.
            h *= max_scale_factor.min(
                min_scale_factor
                    .max(scale_factor * (1. / error_norm).powf(1. / (self.table.p) as f64)),
            );
        }

        self.h = h;
        self.y = self._t.y_next.clone();
        self.t = t_next;

        Some((self.t, self.y.clone()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{Array, Array1};
    use std::fs::File;
    use std::io::Write;

    fn brusselator(_t: &f64, y: &Array1<f64>) -> Array1<f64> {
        Array::from_vec(vec![
            1. + y[0].powi(2) * y[1] - 4. * y[0],
            3. * y[0] - y[0].powi(2) * y[1],
        ])
    }

    #[test]
    fn test_brusselator_adaptive() {
        let t0: f64 = 0.;
        let y0 = Array::from_vec(vec![1.5, 3.]);
        let absolute_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);
        let relative_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);

        let mut integrator = RungeKutta::new_fehlberg(
            t0,
            y0,
            brusselator,
            |t, _| t > &40.,
            0.025,
            relative_tol,
            absolute_tol,
            true,
        );

        integrator.set_step_size(&integrator.guess_initial_step());

        let mut buffer = File::create("examples/brusselator.csv").unwrap();

        for item in integrator {
            writeln!(
                &mut buffer,
                "{:.4}, {:.4}, {:.4}",
                item.0, item.1[0], item.1[1]
            )
            .unwrap();
        }
    }

    #[test]
    fn test_brusselator_fixed() {
        let t0: f64 = 0.;
        let y0 = Array::from_vec(vec![1.5, 3.]);
        let absolute_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);
        let relative_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);

        let integrator = RungeKutta::new_ralston(
            t0,
            y0,
            brusselator,
            |t, _| t > &20.,
            0.025,
            relative_tol,
            absolute_tol,
            false,
        );

        let mut buffer = File::create("examples/brusselator.csv").unwrap();

        for item in integrator {
            writeln!(
                &mut buffer,
                "{:.4}, {:.4}, {:.4}",
                item.0, item.1[0], item.1[1]
            )
            .unwrap();
        }
    }
}
