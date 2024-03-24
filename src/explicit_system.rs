use crate::aux;
use crate::tables;
use ndarray::Array1;
use paste::paste;

macro_rules! impl_new_rksystem {
    ($name:ident) => {
        paste! {
            #[doc="Instantiate an ODE system integrator which uses the " $name:camel " method."]
            pub fn [<new_ $name:lower >] (
                t: f64,
                y: Array1<f64>,
                f: F, predicate: P,
                h: f64
            ) -> Self {
                RungeKuttaSystem {
                    t, y: y.clone(), f, predicate, h,
                    table: tables::[<$name:upper _BT>],
                    _num_variables: y.len(),
                }
            }
        }
    };
}

macro_rules! impl_new_rk_system_adaptive {
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
                absolute_tol: Array1<f64>
            ) -> Self {
                RungeKuttaSystemAdaptive {
                    t, y: y.clone(), f, predicate, h, relative_tol, absolute_tol,
                    table: tables::[<$name:upper _BT_ADAPTIVE>],
                    _num_variables: y.len(),
                }
            }
        }
    }
}

/// Struct for constant step size Runge-Kutta integration. Stores the integration
/// state at every iteration. Implements [`Iterator`].
pub struct RungeKuttaSystem<'a, F, P>
where
    F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
    P: Fn(&f64, &Array1<f64>) -> bool,
{
    t: f64,
    y: Array1<f64>,
    f: F,
    predicate: P,
    h: f64,
    table: tables::ButcherTableau<'a>,
    _num_variables: usize,
}

impl<F, P> RungeKuttaSystem<'_, F, P>
where
    F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
    P: Fn(&f64, &Array1<f64>) -> bool,
{
    impl_new_rksystem!(euler);
    impl_new_rksystem!(ralston);
}

impl<F, P> Iterator for RungeKuttaSystem<'_, F, P>
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
        let mut k: Vec<Array1<f64>> = vec![(self.f)(&self.t, &self.y); self.table.s];

        #[allow(unused_assignments)]
        let mut t: f64 = 0.;

        #[allow(unused_assignments)]
        let mut y = Array1::<f64>::zeros(self._num_variables);

        #[allow(unused_assignments)]
        let mut indx: usize = 0;

        let mut sum = Array1::<f64>::zeros(self._num_variables);

        for i in 1..self.table.s {
            indx = i * (i - 1) / 2;
            t = self.t + self.table.c[i] * self.h;
            y = &self.y + self.h * aux::sum_product_array(&self.table.a[indx..indx + 1], &k[..i]);
            k[i] = (self.f)(&t, &y);
            sum = sum + self.table.b[i] * &k[i];
        }

        self.y = &self.y + self.h * self.table.b[0] * &k[0] + self.h * sum;
        self.t += self.h;

        Some((self.t, self.y.clone()))
    }
}

/// Struct for adaptive step size Runge-Kutta integration. Stores the integration state at every
/// iteration. Implements [`Iterator`].
pub struct RungeKuttaSystemAdaptive<'a, F, P>
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
    table: tables::ButcherTableauAdaptive<'a>,
    _num_variables: usize,
}

///
impl<F, P> RungeKuttaSystemAdaptive<'_, F, P>
where
    F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
    P: Fn(&f64, &Array1<f64>) -> bool,
{
    impl_new_rk_system_adaptive!(Fehlberg);
    impl_new_rk_system_adaptive!(Hueneuler);
    impl_new_rk_system_adaptive!(DormandPrince);

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

impl<F, P> Iterator for RungeKuttaSystemAdaptive<'_, F, P>
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
        let mut k: Vec<Array1<f64>> = vec![(self.f)(&self.t, &self.y); self.table.s];

        #[allow(unused_assignments)]
        let mut t: f64 = 0.;

        #[allow(unused_assignments)]
        let mut y = Array1::<f64>::zeros(self._num_variables);

        #[allow(unused_assignments)]
        let mut indx: usize = 0;

        let mut h = self.h;
        let mut error_norm = 2.;
        let mut t_next: f64 = 0.;
        let mut y_next = Array1::<f64>::zeros(self._num_variables);

        #[allow(unused_assignments)]
        let mut y_lower = Array1::<f64>::zeros(self._num_variables);

        let mut sum = Array1::<f64>::zeros(self._num_variables);

        let scale_factor: f64 = 0.9;
        let min_scale_factor: f64 = 0.1;
        let max_scale_factor: f64 = 2.;

        let mut count: usize = 0;
        while error_norm > 1. {
            count += 1;

            sum.fill(0.);

            for i in 1..self.table.s {
                indx = i * (i - 1) / 2;
                t = self.t + self.table.c[i] * self.h;
                y = &self.y
                    + self.h * aux::sum_product_array(&self.table.a[indx..indx + 1], &k[..i]);
                k[i] = (self.f)(&t, &y);
                sum = sum + self.table.b[i] * &k[i];
            }

            t_next = self.t + h;
            y_next = &self.y + h * self.table.b[0] * &k[0] + h * &sum;
            y_lower = &self.y + h * aux::sum_product_array(&self.table.b2, &k);

            error_norm = self.calc_error_norm(&y_next, &y_lower);

            h *= max_scale_factor.min(
                min_scale_factor
                    .max(scale_factor * (1. / error_norm).powf(1. / (self.table.p) as f64)),
            );

            // println!(
            //     "count = {:}, h = {:?}, t = {:.4}, y = {:.4}, norm = {:.4}, scale = {:.4}",
            //     count, h, t_next, y_next, error_norm, _tmp
            // );
        }
        // println!("---- Advancing to next");
        self.h = h;
        self.y = y_next;
        self.t = t_next;

        Some((self.t, self.y.clone()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{s, Array, Array1};
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn test_euler_system() {
        let y0 = Array1::<f64>::zeros(3);

        let mut integrator = RungeKuttaSystem::new_euler(
            0.,
            y0,
            |_, _| Array1::<f64>::from_elem(3, 2.),
            |t, _| t > &5.,
            1.0,
        );

        for item in integrator {
            println!("{:?}", item);
        }
    }

    #[test]
    fn test_ralston_system() {
        let y0 = Array1::<f64>::zeros(3) + 1.;

        let mut integrator = RungeKuttaSystem::new_ralston(
            1.,
            y0,
            |_, y| y.map(|x| x.tan()) + 1.,
            |t, _| t > &1.075,
            0.025,
        );

        for item in integrator {
            println!("{:?}", item);
        }
    }

    #[test]
    fn test_fehlberg_system() {
        fn brusselator(t: &f64, y: &Array1<f64>) -> Array1<f64> {
            Array::from_vec(vec![
                1. + y[0].powi(2) * y[1] - 4. * y[0],
                3. * y[0] - y[0].powi(2) * y[1],
            ])
        }

        let t0: f64 = 0.;
        let y0 = Array::from_vec(vec![1.5, 3.]);
        let absolute_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);
        let relative_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);

        let mut integrator = RungeKuttaSystemAdaptive::new_fehlberg(
            t0,
            y0,
            brusselator,
            |t, _| t > &20.,
            0.025,
            relative_tol,
            absolute_tol,
        );

        integrator.h = integrator.guess_initial_step();
        // println!("{:?}", integrator.h);
        // panic!();

        for (i, item) in integrator.enumerate() {
            println!("{:?}", item);
            if i > 5 {
                break;
            }
        }
    }

    #[test]
    fn test_brusselator_adaptive() {
        fn brusselator(t: &f64, y: &Array1<f64>) -> Array1<f64> {
            Array::from_vec(vec![
                1. + y[0].powi(2) * y[1] - 4. * y[0],
                3. * y[0] - y[0].powi(2) * y[1],
            ])
        }

        let t0: f64 = 0.;
        let y0 = Array::from_vec(vec![1.5, 3.]);
        let absolute_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);
        let relative_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);

        let mut integrator = RungeKuttaSystemAdaptive::new_fehlberg(
            t0,
            y0,
            brusselator,
            |t, _| t > &20.,
            0.025,
            relative_tol,
            absolute_tol,
        );

        integrator.h = integrator.guess_initial_step();

        let mut buffer = File::create("brusselator.csv").unwrap();

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
        fn brusselator(t: &f64, y: &Array1<f64>) -> Array1<f64> {
            Array::from_vec(vec![
                1. + y[0].powi(2) * y[1] - 4. * y[0],
                3. * y[0] - y[0].powi(2) * y[1],
            ])
        }

        let t0: f64 = 0.;
        let y0 = Array::from_vec(vec![1.5, 3.]);
        let absolute_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);
        let relative_tol = Array::from_vec(vec![1.0e-4, 1.0e-4]);

        let mut integrator =
            RungeKuttaSystem::new_ralston(t0, y0, brusselator, |t, _| t > &20., 0.025);

        let mut buffer = File::create("brusselator.csv").unwrap();

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
