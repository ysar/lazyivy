use crate::tables;
use ndarray::{array, s, Array1, Array2};

/// Builder struct for `RungeKutta`. It can be called directly or via `RungeKutta::builder()`. Like
/// `RungeKutta`, it is generic over the evaluation function and predicate.
pub struct RungeKuttaBuilder<F, P>
where
    F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
    P: Fn(&f64, &Array1<f64>) -> bool,
{
    t: f64,
    y: Array1<f64>,
    f: F,
    predicate: P,
    step: f64,
    relative_tol: Array1<f64>,
    absolute_tol: Array1<f64>,
    do_adaptive: bool,
    max_step_size: f64,
    method: String,
}

impl<F, P> RungeKuttaBuilder<F, P>
where
    F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
    P: Fn(&f64, &Array1<f64>) -> bool,
{
    /// Create a new `RungeKuttaBuilder` instance.
    pub fn new(f: F, predicate: P) -> Self {
        // Create a default RungeKuttaBuilder. Cannot do this via the Default trait because of
        // generics.
        RungeKuttaBuilder {
            t: 0.,
            y: array![0.],
            f,
            predicate,
            step: 1.,
            relative_tol: array![1.0e-4],
            absolute_tol: array![1.0e-4],
            do_adaptive: false,
            max_step_size: f64::INFINITY,
            method: "Euler".to_string(),
        }
    }

    /// Set the initial condition `(t0, y0)` (Default: `t0 = 0.` and `y0 = array1![0.]`)
    pub fn initial_condition(mut self, t: f64, y: Array1<f64>) -> Self {
        self.t = t;
        self.y = y;
        self
    }

    /// Sets the initial step size of the iteration.
    pub fn initial_step_size(mut self, step: f64) -> Self {
        self.step = step;
        self
    }

    /// Sets the Runge-Kutta method to be used. Options are :
    /// `euler`, `ralston`, `hueneuler`, `bogackishampine`, `fehlberg`, `dormandprince`.
    pub fn method(mut self, method: &str, do_adaptive: bool) -> Self {
        self.method = method.to_string();
        self.do_adaptive = do_adaptive;
        self
    }

    /// Sets the relative and absolute tolerances for the iterations. Both must be arrays of the
    /// same size as the vector of unknowns and the result of the evaluation function.
    pub fn tolerances(mut self, absolute_tol: Array1<f64>, relative_tol: Array1<f64>) -> Self {
        self.relative_tol = relative_tol;
        self.absolute_tol = absolute_tol;
        self
    }

    /// Species whether to use adaptive step size modification based on a lower order method. Not
    /// enabled by default.
    pub fn do_adaptive(mut self, is_adaptive: bool) -> Self {
        self.do_adaptive = is_adaptive;
        self
    }

    /// Specifies the maximum step-size for adaptive step-size problems. This is sometimes very
    /// useful if the tolerances are not well known beforehand.
    pub fn set_max_step_size(mut self, max_step: f64) -> Self {
        self.max_step_size = max_step;
        self
    }

    /// Consumes the `RungeKuttaBuilder` and returns a `Result` containing the `RungeKutta` struct.
    /// This needs to be called at the very end of the construction chain.
    pub fn build(self) -> Result<RungeKutta<F, P>, String> {
        let table = tables::get_rungekutta_coefficients(&self.method)?;

        let num_variables: usize = self.y.len();
        let f0 = (self.f)(&self.t, &self.y);

        if f0.len() != num_variables {
            return Err("Evaluation function inconsistent with initial condition.".to_string());
        }

        if self.relative_tol.len() != num_variables || self.absolute_tol.len() != num_variables {
            return Err(
                "Tolerances need to be arrays with the same size as the vector of unknowns."
                    .to_string(),
            );
        }

        match self.method.as_str() {
            "euler" => {
                return Err(format!(
                    "Adaptive step-size not implemented for method : '{:}'",
                    self.method
                ))
            }
            "ralston" => {
                return Err(format!(
                    "Adaptive step-size not implemented for method : '{:}'",
                    self.method
                ))
            }
            _ => {}
        }

        // Pre-allocate some arrays used in the iterations. Cloning here is fine, we are only doing
        // it once at the beginning.

        let tmp_variables = _TempArrays::new(num_variables, table.num_stages);

        Ok(RungeKutta {
            t: self.t,
            y: self.y,
            f: self.f,
            predicate: self.predicate,
            step: self.step,
            relative_tol: self.relative_tol,
            absolute_tol: self.absolute_tol,
            do_adaptive: self.do_adaptive,
            table,
            max_step_size: self.max_step_size,
            _num_variables: num_variables,
            _t: tmp_variables,
        })
    }
}

/// Struct to store preallocated arrays
struct _TempArrays {
    k: Array2<f64>,
    y: Array1<f64>,
    y_next: Array1<f64>,
    y_err: Array1<f64>,
    sum: Array1<f64>,
}

impl _TempArrays {
    fn new(num_variables: usize, num_stages: usize) -> Self {
        let zero1 = Array1::<f64>::zeros(num_variables);
        let zero2 = Array2::<f64>::zeros((num_stages, num_variables));

        _TempArrays {
            k: zero2.clone(),
            y: zero1.clone(),
            y_next: zero1.clone(),
            y_err: zero1.clone(),
            sum: zero1.clone(),
        }
    }
}

/// Struct for adaptive step size Runge-Kutta integration. Stores the integration state at every
/// iteration. Implements [`Iterator`].
pub struct RungeKutta<F, P>
where
    F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
    P: Fn(&f64, &Array1<f64>) -> bool,
{
    t: f64,
    y: Array1<f64>,
    f: F,
    predicate: P,
    step: f64,
    relative_tol: Array1<f64>,
    absolute_tol: Array1<f64>,
    do_adaptive: bool,
    table: tables::ButcherTableau,
    max_step_size: f64,
    _num_variables: usize,
    _t: _TempArrays,
}

///
impl<F, P> RungeKutta<F, P>
where
    F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
    P: Fn(&f64, &Array1<f64>) -> bool,
{
    /// Returns a `RungeKuttaBuilder` for building a `RungeKutta` struct.
    pub fn builder(f: F, predicate: P) -> RungeKuttaBuilder<F, P> {
        RungeKuttaBuilder::new(f, predicate)
    }

    /// Set step-size.
    pub fn set_step_size(&mut self, step: &f64) {
        self.step = *step;
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
            (0.01 / d1.max(d2)).powf(1. / (self.table.order + 1) as f64)
        };

        h1.min(100. * h0)
    }
}

impl<F, P> Iterator for RungeKutta<F, P>
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

        self._t
            .k
            .rows_mut()
            .into_iter()
            .for_each(|mut row| row.assign(&f_now));
        self._t.y.fill(0.);
        self._t.y_next.fill(0.);
        self._t.y_err.fill(0.);

        #[allow(unused_assignments)]
        let mut t: f64 = 0.;

        #[allow(unused_assignments)]
        let mut indx: usize = 0;

        #[allow(unused_assignments)]
        let mut t_next: f64 = 0.;

        let mut step = self.step;

        // Set the error norm to be larger than 1. before first while condition check.
        let mut error_norm = 2.;

        // Scale factors for changing the step-size. Prevents overzealous changes.
        let scale_factor: f64 = 0.9;
        let min_scale_factor: f64 = 0.1;
        let max_scale_factor: f64 = 2.;

        while error_norm > 1. {
            self._t.sum.fill(0.);

            for i_stage in 1..self.table.num_stages {
                // Index into the lower triangular matrix.
                indx = i_stage * (i_stage - 1) / 2;

                // Advance time in the inner stage.
                t = self.t + self.table.c[i_stage] * self.step;

                // Calculate `y` in the inner stage using k[0], k[1], ... k[stage - 1].
                self._t.y = &self.y
                    + self.step
                        * self
                            .table
                            .a
                            .slice(s![indx..indx + i_stage])
                            .dot(&self._t.k.slice(s![..i_stage, ..]));

                // Calculate `k` for the next stage.
                self._t.k.row_mut(i_stage).assign(&(self.f)(&t, &self._t.y));

                // Add the stage result and move to next stage.
                self._t.sum = &self._t.sum + self.table.b[i_stage] * &self._t.k.row(i_stage);
            }

            // Advance time and output in the outer stage.
            t_next = self.t + step;
            self._t.y_next =
                &self.y + step * self.table.b[0] * &self._t.k.row(0) + step * &self._t.sum;

            // If no adaptive step-size requested, we can break out of the loop now.
            if !self.do_adaptive {
                break;
            }

            // Calculate the result using the error estimator coefficients.
            self._t.y_err = &self.y + step * self.table.b2.dot(&self._t.k);

            error_norm = self.calc_error_norm(&self._t.y_next, &self._t.y_err);

            // Scale step-size based on the norm of the error between the result and the error
            // estimator.
            step *=
                max_scale_factor
                    .min(min_scale_factor.max(
                        scale_factor * (1. / error_norm).powf(1. / (self.table.order) as f64),
                    ));

            if step > self.max_step_size {
                step = self.max_step_size;
            }
        }

        self.step = step;
        self.y = self._t.y_next.clone();
        self.t = t_next;

        Some((self.t, self.y.clone()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{array, Array1};

    fn brusselator(_t: &f64, y: &Array1<f64>) -> Array1<f64> {
        array![
            1. + y[0].powi(2) * y[1] - 4. * y[0],
            3. * y[0] - y[0].powi(2) * y[1],
        ]
    }

    #[test]
    fn test_brusselator_fixed() -> Result<(), String> {
        let t0: f64 = 0.;
        let y0 = array![1.5, 3.];
        let absolute_tol = array![1.0e-4, 1.0e-4];
        let relative_tol = array![1.0e-4, 1.0e-4];

        let integrator = RungeKutta::builder(brusselator, |t, _| *t > 40.)
            .initial_condition(t0, y0)
            .initial_step_size(0.025)
            .method("fehlberg", true)
            .tolerances(absolute_tol, relative_tol)
            .set_max_step_size(0.25)
            .build()?;

        for _item in integrator {}
        Ok(())
    }
}
