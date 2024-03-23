use crate::aux;
use crate::tables;
use ndarray::{Array1};

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
}

impl<F, P> RungeKuttaSystem<'_, F, P>
where
    F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
    P: Fn(&f64, &Array1<f64>) -> bool,
{
    /// Instantiate an integrator which uses the Euler method for a system 
    /// of ODEs.
    pub fn new_euler(t0: f64, y0: Array1<f64>, f_in: F, p: P, h_in: f64) -> Self {
        RungeKuttaSystem {
            t: t0,
            y: y0,
            f: f_in,
            predicate: p,
            h: h_in,
            table: tables::EULER_BT,
        }
    }

    /// Instantiate an integrator which uses the Euler method for a system 
    /// of ODEs.
    pub fn new_ralston(t0: f64, y0: Array1<f64>, f_in: F, fstop: P, h_in: f64) -> Self {
        RungeKuttaSystem {
            t: t0,
            y: y0,
            f: f_in,
            predicate: fstop,
            h: h_in,
            table: tables::RALSTON_BT,
        }
    }
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

        let mut t: f64 = 0.;
        let mut y = Array1::<f64>::zeros(self.y.len());

        // See butcher.rs for note on indexing logic for Butcher table

        let mut indx: usize = 0;

        let mut sum = Array1::<f64>::zeros(self.y.len());

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

// pub struct RungeKuttaAdaptive<'a, F, P>
// where
//     F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
//     P: Fn(&f64, &Array1<f64>) -> bool,
// {
//     t: f64,
//     y: f64,
//     f: F,
//     predicate: P,
//     h: f64,
//     err0: f64,
//     table: tables::ButcherTableauAdaptive<'a>,
// }

// impl<F, P> RungeKuttaAdaptive<'_, F, P>
// where
//     F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
//     P: Fn(&f64, &Array1<f64>) -> bool,
// {
//     pub fn new_fehlberg(t0: f64, y0: f64, f_in: F, fstop: P, h_in: f64, err_in: f64) -> Self {
//         RungeKuttaAdaptive {
//             t: t0,
//             y: y0,
//             f: f_in,
//             predicate: fstop,
//             h: h_in,
//             err0: err_in,
//             table: tables::FEHLBERG_BT_ADAPTIVE,
//         }
//     }
// }

// impl<F, P> RungeKuttaAdaptive<'_, F, P>
// where
//     F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
//     P: Fn(&f64, &Array1<f64>) -> bool,
// {
//     pub fn guess_initial_step(&self) -> f64 {
//         let f0 = (self.f)(&self.t, &self.y);
//         let d0 = (self.y / self.err0).abs();
//         let d1 = (f0 / self.err0).abs();

//         let h0: f64 = if d0 < 1.0e-5 || d1 < 1.0e-5 {
//             1.0e-6
//         } else {
//             0.01 * (d0 / d1)
//         };

//         let y1 = self.y + h0 * f0;
//         let f1 = (self.f)(&(self.t + h0), &y1);

//         let d2 = ((f1 - f0) / self.err0).abs() / h0;

//         let h1: f64 = if d1.max(d2) < 1e-15 {
//             (h0 * 1.0e-3).max(1.0e-6)
//         } else {
//             (0.01 / d1.max(d2)).powf(1. / (self.table.p + 1) as f64)
//         };

//         h1.min(100. * h0)
//     }
// }

// impl<F, P> Iterator for RungeKuttaAdaptive<'_, F, P>
// where
//     F: Fn(&f64, &Array1<f64>) -> Array1<f64>,
//     P: Fn(&f64, &Array1<f64>) -> bool,
// {
//     type Item = (f64, f64);

//     fn next(&mut self) -> Option<Self::Item> {
//         // Check stop condition at every iteration
//         if (self.predicate)(&self.t, &self.y) {
//             return None;
//         }

//         // RK Logic

//         // Store k[i] values in a vector, initialize with the first value
//         let mut k: Vec<f64> = vec![(self.f)(&self.t, &self.y); self.table.s];

//         let mut t: f64 = 0.;
//         let mut y: f64 = 0.;
//         let mut indx: usize = 0;
//         let mut h = self.h;
//         let mut err = self.err0 + 1.;
//         let mut y_next: f64 = 0.;
//         let mut y_err: f64;
//         let mut t_next: f64 = 0.;

//         let facmax: f64 = 2.;
//         let fac: f64 = 0.9;
//         let facmin: f64 = 0.;

//         while err > self.err0 {
//             println!("h = {:.3}, t = {:.3}, y = {:.3}", h, t_next, y_next);

//             y_next = self.y
//                 + h * self.table.b[0] * k[0]
//                 + h * (1..self.table.s)
//                     .map(|i| {
//                         indx = i * (i - 1) / 2;

//                         t = self.t + self.table.c[i] * h;
//                         y = self.y + h * aux::sum_product(&self.table.a[indx..indx + i], &k[..i]);

//                         k[i] = (self.f)(&t, &y);

//                         self.table.b[i] * k[i]
//                     })
//                     .sum::<f64>();

//             t_next = self.t + h;

//             y_err = self.y + h * aux::sum_product(&self.table.b2, &k);

//             err = (y_next - y_err).abs();

//             h *= facmax
//                 .min(facmin.max(fac * (self.err0 / err).powf(1. / (self.table.p + 1) as f64)));
//             // println!("h = {:.3}, t = {:.3}, y = {:.3}", h, t_next, y_next);
//         }
//         // println!("---- Advancing to next");
//         self.h = h;
//         self.y = y_next;
//         self.t = t_next;

//         Some((self.t, self.y))
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array1;

    #[test]
    fn test_euler_system() {

        let y0 = Array1::<f64>::zeros(3);

        let mut integrator = RungeKuttaSystem::new_euler(
            0.,
            y0, 
            |_, _| Array1::<f64>::from_elem(3, 2.), 
            |t, _| t > &5., 
            1.0);
        
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
            |_, y| { y.map(|x| x.tan()) + 1.} , 
            |t, _| t > &1.075, 
            0.025);
        
        for item in integrator {
            println!("{:?}", item);
        }
    }
}