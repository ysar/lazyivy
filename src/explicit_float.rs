use crate::butcher;

struct RungeKuttaIntegrator<'a, F, S>
where
    F: Fn(&f64, &f64) -> f64,
    S: Fn(&f64, &f64) -> bool,
{
    t: f64,
    y: f64,
    f: F,
    stop_now: S,
    h: f64,
    butcher: butcher::ButcherTableau<'a>,
}

impl<F, S> RungeKuttaIntegrator<'_, F, S>
where
    F: Fn(&f64, &f64) -> f64,
    S: Fn(&f64, &f64) -> bool,
{
    pub fn new_euler(t0: f64, y0: f64, f_in: F, fstop: S, h_in: f64) -> Self {
        RungeKuttaIntegrator {
            t: t0,
            y: y0,
            f: f_in,
            stop_now: fstop,
            h: h_in,
            butcher: butcher::EULER_BT,
        }
    }

    pub fn new_ralston(t0: f64, y0: f64, f_in: F, fstop: S, h_in: f64) -> Self {
        RungeKuttaIntegrator {
            t: t0,
            y: y0,
            f: f_in,
            stop_now: fstop,
            h: h_in,
            butcher: butcher::RALSTON_BT,
        }
    }
}

impl<F, S> Iterator for RungeKuttaIntegrator<'_, F, S>
where
    F: Fn(&f64, &f64) -> f64,
    S: Fn(&f64, &f64) -> bool,
{
    type Item = (f64, f64);

    fn next(&mut self) -> Option<Self::Item> {
        // Check stop condition at every iteration
        if (self.stop_now)(&self.t, &self.y) {
            return None;
        }

        // RK Logic

        // Store k[i] values in a vector, initialize with the first value
        let mut k: Vec<f64> = vec![(self.f)(&self.t, &self.y); self.butcher.s];

        let mut t: f64 = 0.;
        let mut y: f64 = 0.;

        // See butcher.rs for note on indexing logic for Butcher table
        
        let mut indx: usize = 0;

        self.t += self.h;
        self.y += 
              self.h * self.butcher.b[0] * k[0] 
            + self.h * (1..self.butcher.s).map(|i| 
                {
                    indx = i * (i - 1) / 2; 

                    t = self.t + self.butcher.c[i] * self.h;
                    
                    y = self.y + self.h 
                        * (0..i)
                            .map(|j| self.butcher.a[indx + j] * k[j])
                            .sum::<f64>();
                    
                    k[i] = (self.f)(&t, &y);
                    
                    self.butcher.b[i] * k[i]
                }
            ).sum::<f64>();

        Some((self.t, self.y))
    }
}


struct RungeKuttaIntegratorAdaptive<'a, F, S>
where
    F: Fn(&f64, &f64) -> f64,
    S: Fn(&f64, &f64) -> bool,
{
    t: f64,
    y: f64,
    f: F,
    stop_now: S,
    h: f64,
    butcher: butcher::ButcherTableauAdaptive<'a>,
}

impl<F, S> RungeKuttaIntegratorAdaptive<'_, F, S>
where
    F: Fn(&f64, &f64) -> f64,
    S: Fn(&f64, &f64) -> bool,
{
    pub fn new_fehlberg(t0: f64, y0: f64, f_in: F, fstop: S, h_in: f64) -> Self {
        RungeKuttaIntegratorAdaptive {
            t: t0,
            y: y0,
            f: f_in,
            stop_now: fstop,
            h: h_in,
            butcher: butcher::FEHLBERG_BT_ADAPTIVE,
        }
    }
}

impl<F, S> Iterator for RungeKuttaIntegratorAdaptive<'_, F, S>
where
    F: Fn(&f64, &f64) -> f64,
    S: Fn(&f64, &f64) -> bool,
{
    type Item = (f64, f64);

    fn next(&mut self) -> Option<Self::Item> {
        // Check stop condition at every iteration
        if (self.stop_now)(&self.t, &self.y) {
            return None;
        }

        // RK Logic

        // Store k[i] values in a vector, initialize with the first value
        let mut k: Vec<f64> = vec![(self.f)(&self.t, &self.y); self.butcher.s + 1];

        let mut t: f64 = 0.;
        let mut y: f64 = 0.;
        

        // Note: 
        // The indexing logic into the a field of a Butcher table (explicit) 
        // stored in lower triangular form is 
        // indx = i * (i - 1) / 2 + j


        let mut indx: usize = 0;

        self.t += self.h;
        self.y += 
              self.h * self.butcher.b[0] * k[0] 
            + self.h * (1..self.butcher.s).map(|i| 
                {
                    indx = i * (i - 1) / 2; 

                    t = self.t + self.butcher.c[i] * self.h;
                    
                    y = self.y + self.h 
                        * (0..i)
                            .map(|j| self.butcher.a[indx + j] * k[j])
                            .sum::<f64>();
                    
                    k[i] = (self.f)(&t, &y);
                    
                    self.butcher.b[i] * k[i]
                }
            ).sum::<f64>();

        Some((self.t, self.y))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_euler_constant() {
        let mut integrator = RungeKuttaIntegrator::new_euler(
            0., 0., |_, _| 2., |_, y| y > &5., 1.0);

        let result_correct = vec![(1.0, 2.0), (2.0, 4.0), (3.0, 6.0)];
        assert_eq!(result_correct, integrator.collect::<Vec<_>>());
    }

    #[test]
    fn test_ralston() {
        let integrator = RungeKuttaIntegrator::new_ralston(
            1.,
            1.,
            |_, y| y.tan() + 1.,
            |t, _| t > &1.075,
            0.025,
        );

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
        let integrator = RungeKuttaIntegratorAdaptive::new_fehlberg(
            1.,
            1.,
            |_, y| y.tan() + 1.,
            |t, _| t > &1.075,
            0.025,
        );

        let result_correct = vec![
            (1.025, 1.06697),
            (1.05, 1.141636),
            (1.075, 1.22822),
            (1.1, 1.33786),
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
}