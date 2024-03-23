use crate::aux;
use crate::tables;

/// Struct for constant step size Runge-Kutta integration. Stores the integration
/// state at every iteration. Implements [`Iterator`].
///
/// # Usage:
/// Instantiate an [`RungeKutta`] instance using coefficients
/// for the Euler method. `t0` and `y0` are `f64` and represent the
/// initial condition. Two function pointers `F` and `P` need to be passed by
/// by the user. `F` is the evaluation function, i.e., the right hand side
/// function for an ODE `dy/dt = f(t, y)`. The second function pointer `P` is
/// the predicate that dictates the stop condition.
/// `F` and `P`function pointers that implement the [`Fn`] trait, so they
///  can also be closures that do not capture their environment, s
/// since these closures can be coerced into function pointers.
/// ```
/// // Create integrator to solve dy/dt=2 starting from t=0, y=0 until
/// // y exceeds y=5, using constant step size of h=1.
/// # use lazyivy::explicit_single::RungeKutta;
/// let integrate = RungeKutta::new_euler(
///     0., 0., |_, _| 2., |_, y| y > &5., 1.0);
/// ```
/// This will create an iterator but will not consume it since iterators
/// are lazy. To consume the iterator, you have various choices.  
/// You can collect all integration steps into a vector.
/// ```
/// # use lazyivy::explicit_single::RungeKutta;
/// # let integrate = RungeKutta::new_euler(
/// #     0., 0., |_, _| 2., |_, y| y > &5., 1.0);
/// let result_all = integrate.collect::<Vec<_>>();
///   // This creates a vector of tuples vec![(t0, y0), (t1, y1) ... (tN, yN)]
/// ```
/// Or you can iterate till the last value and only keep that.
/// ```
/// # use lazyivy::explicit_single::RungeKutta;
/// # let integrate = RungeKutta::new_euler(
/// #    0., 0., |_, _| 2., |_, y| y > &5., 1.0);
/// let result_last = integrate.last();
///   // result_last will be (tN, yN)
/// ```
/// Or you can use any of the methods implemented within the [`Iterator`] trait
/// e.g. `map`, `for_each`, etc.  
/// Various Runge-Kutta methods with varying number of stages and accuracy are
/// provided and can be used similar to the example after instantiating with a
/// call to `RungeKutta::new_{name-of-rk-method}`.
pub struct RungeKutta<'a, F, P>
where
    F: Fn(&f64, &f64) -> f64,
    P: Fn(&f64, &f64) -> bool,
{
    t: f64,
    y: f64,
    f: F,
    predicate: P,
    h: f64,
    table: tables::ButcherTableau<'a>,
}

impl<F, P> RungeKutta<'_, F, P>
where
    F: Fn(&f64, &f64) -> f64,
    P: Fn(&f64, &f64) -> bool,
{
    /// Instantiate an integration using the Euler method
    pub fn new_euler(t0: f64, y0: f64, f_in: F, p: P, h_in: f64) -> Self {
        RungeKutta {
            t: t0,
            y: y0,
            f: f_in,
            predicate: p,
            h: h_in,
            table: tables::EULER_BT,
        }
    }

    /// Instantiate an integration using the Ralston method
    pub fn new_ralston(t0: f64, y0: f64, f_in: F, fstop: P, h_in: f64) -> Self {
        RungeKutta {
            t: t0,
            y: y0,
            f: f_in,
            predicate: fstop,
            h: h_in,
            table: tables::RALSTON_BT,
        }
    }
}

impl<F, P> Iterator for RungeKutta<'_, F, P>
where
    F: Fn(&f64, &f64) -> f64,
    P: Fn(&f64, &f64) -> bool,
{
    type Item = (f64, f64);

    fn next(&mut self) -> Option<Self::Item> {
        // Check stop condition at every iteration
        if (self.predicate)(&self.t, &self.y) {
            return None;
        }

        // RK Logic

        // Store k[i] values in a vector, initialize with the first value
        let mut k: Vec<f64> = vec![(self.f)(&self.t, &self.y); self.table.s];

        let mut t: f64 = 0.;
        let mut y: f64 = 0.;

        // See butcher.rs for note on indexing logic for Butcher table

        let mut indx: usize = 0;

        self.y += self.h * self.table.b[0] * k[0]
            + self.h
                * (1..self.table.s)
                    .map(|i| {
                        indx = i * (i - 1) / 2;

                        t = self.t + self.table.c[i] * self.h;
                        y = self.y
                            + self.h * aux::sum_product(&self.table.a[indx..indx + i], &k[..i]);

                        k[i] = (self.f)(&t, &y);

                        self.table.b[i] * k[i]
                    })
                    .sum::<f64>();

        self.t += self.h;

        Some((self.t, self.y))
    }
}

/// Struct for adaptive step-size Runge-Kutta integration. Makes use of a lower
/// order accurate scheme to calculate the error and scales the step size 
/// accordingly.
pub struct RungeKuttaAdaptive<'a, F, P>
where
    F: Fn(&f64, &f64) -> f64,
    P: Fn(&f64, &f64) -> bool,
{
    t: f64,
    y: f64,
    f: F,
    predicate: P,
    h: f64,
    err0: f64,
    table: tables::ButcherTableauAdaptive<'a>,
}

impl<F, P> RungeKuttaAdaptive<'_, F, P>
where
    F: Fn(&f64, &f64) -> f64,
    P: Fn(&f64, &f64) -> bool,
{   
    /// Instantiates an integration using the adaptive Fehlberg method
    pub fn new_fehlberg(t0: f64, y0: f64, f_in: F, fstop: P, h_in: f64, err_in: f64) -> Self {
        RungeKuttaAdaptive {
            t: t0,
            y: y0,
            f: f_in,
            predicate: fstop,
            h: h_in,
            err0: err_in,
            table: tables::FEHLBERG_BT_ADAPTIVE,
        }
    }
}

impl<F, P> RungeKuttaAdaptive<'_, F, P>
where
    F: Fn(&f64, &f64) -> f64,
    P: Fn(&f64, &f64) -> bool,
{   
    /// Provides an initial guess for adaptive step size algorithms. Follows
    /// the algorithm written in the book by Harrier, Nørsett, Wanner.
    pub fn guess_initial_step(&self) -> f64 {
        let f0 = (self.f)(&self.t, &self.y);
        let d0 = (self.y / self.err0).abs();
        let d1 = (f0 / self.err0).abs();

        let h0: f64 = if d0 < 1.0e-5 || d1 < 1.0e-5 {
            1.0e-6
        } else {
            0.01 * (d0 / d1)
        };

        let y1 = self.y + h0 * f0;
        let f1 = (self.f)(&(self.t + h0), &y1);

        let d2 = ((f1 - f0) / self.err0).abs() / h0;

        let h1: f64 = if d1.max(d2) < 1e-15 {
            (h0 * 1.0e-3).max(1.0e-6)
        } else {
            (0.01 / d1.max(d2)).powf(1. / (self.table.p + 1) as f64)
        };

        h1.min(100. * h0)
    }
}

impl<F, P> Iterator for RungeKuttaAdaptive<'_, F, P>
where
    F: Fn(&f64, &f64) -> f64,
    P: Fn(&f64, &f64) -> bool,
{
    type Item = (f64, f64);

    fn next(&mut self) -> Option<Self::Item> {
        // Check stop condition at every iteration
        if (self.predicate)(&self.t, &self.y) {
            return None;
        }

        // RK Logic

        // Store k[i] values in a vector, initialize with the first value
        let mut k: Vec<f64> = vec![(self.f)(&self.t, &self.y); self.table.s];

        let mut t: f64 = 0.;
        let mut y: f64 = 0.;
        let mut indx: usize = 0;
        let mut h = self.h;
        let mut err = self.err0 + 1.;
        let mut y_next: f64 = 0.;
        let mut y_err: f64;
        let mut t_next: f64 = 0.;

        let facmax: f64 = 2.;
        let fac: f64 = 0.9;
        let facmin: f64 = 0.;

        while err > self.err0 {
            println!("h = {:.3}, t = {:.3}, y = {:.3}", h, t_next, y_next);

            y_next = self.y
                + h * self.table.b[0] * k[0]
                + h * (1..self.table.s)
                    .map(|i| {
                        indx = i * (i - 1) / 2;

                        t = self.t + self.table.c[i] * h;
                        y = self.y + h * aux::sum_product(&self.table.a[indx..indx + i], &k[..i]);

                        k[i] = (self.f)(&t, &y);

                        self.table.b[i] * k[i]
                    })
                    .sum::<f64>();

            t_next = self.t + h;

            y_err = self.y + h * aux::sum_product(&self.table.b2, &k);

            err = (y_next - y_err).abs();

            h *= facmax
                .min(facmin.max(fac * (self.err0 / err).powf(1. / (self.table.p + 1) as f64)));
            // println!("h = {:.3}, t = {:.3}, y = {:.3}", h, t_next, y_next);
        }
        // println!("---- Advancing to next");
        self.h = h;
        self.y = y_next;
        self.t = t_next;

        Some((self.t, self.y))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_euler_constant() {
        let mut integrator = RungeKutta::new_euler(0., 0., |_, _| 2., |_, y| y > &5., 1.0);

        let result_correct = vec![(1.0, 2.0), (2.0, 4.0), (3.0, 6.0)];
        assert_eq!(result_correct, integrator.collect::<Vec<_>>());
    }

    #[test]
    fn test_ralston() {
        let integrator =
            RungeKutta::new_ralston(1., 1., |_, y| y.tan() + 1., |t, _| t > &1.075, 0.025);

        let result_correct = vec![
            (1.025, 1.06686),
            (1.05, 1.14133),
            (1.075, 1.22741),
            (1.1, 1.33507),
        ];

        let result = integrator.collect::<Vec<_>>();

        let not_equal_vals = result_correct
            .iter()
            .zip(result.iter())
            .filter(|&(x, y)| (x.1 - y.1).abs() > 0.001)
            .count();

        assert!(
            not_equal_vals == 0,
            "{:?}\n{:?}\n Not equal vals ={:?}",
            result,
            result_correct,
            not_equal_vals
        );
    }

    #[test]
    fn test_fehlberg() {
        let mut integrator =
            RungeKuttaAdaptive::new_fehlberg(1., 1., |t, _| t * t, |t, _| t > &10., 0.025, 0.0001);

        integrator.h = integrator.guess_initial_step();

        for iteration in integrator {
            // println!("{:?}", iteration);
        }
    }
}
