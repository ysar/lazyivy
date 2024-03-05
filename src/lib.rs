struct ButcherTableau<'t> {
    s: usize,
    a: &'t [f64],
    b: &'t [f64],
    c: &'t [f64],
}

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
    butcher: ButcherTableau<'a>,
}

impl<F, S> RungeKuttaIntegrator<'_, F, S>
where
    F: Fn(&f64, &f64) -> f64,
    S: Fn(&f64, &f64) -> bool,
{
    fn new_euler(t0: f64, y0: f64, f_in: F, fstop: S, h_in: f64) -> Self {
        RungeKuttaIntegrator {
            t: t0,
            y: y0,
            f: f_in,
            stop_now: fstop,
            h: h_in,
            butcher: ButcherTableau {
                s: 1,
                a: &[0.],
                b: &[1.],
                c: &[0.],
            },
        }
    }

    fn new_ralston(t0: f64, y0: f64, f_in: F, fstop: S, h_in: f64) -> Self {
        RungeKuttaIntegrator {
            t: t0,
            y: y0,
            f: f_in,
            stop_now: fstop,
            h: h_in,
            butcher: ButcherTableau {
                s: 2,
                a: &[0., 0., 2. / 3., 0.],
                b: &[1. / 4., 3. / 4.],
                c: &[0., 2. / 3.],
            },
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
        let mut k: Vec<f64> = vec![(self.f)(&self.t, &self.y); self.butcher.s + 1];

        let mut t: f64 = 0.;
        let mut y: f64 = 0.;

        self.t += self.h;
        self.y += self.h
            * (0..self.butcher.s)
                .map(|i| {
                    t = self.t + self.butcher.c[i] * self.h;
                    y = self.y
                        + self.h
                            * (0..i)
                                .map(|j| self.butcher.a[i * self.butcher.s + j] * k[j])
                                .sum::<f64>();
                    k[i] = (self.f)(&t, &y);
                    self.butcher.b[i] * k[i]
                })
                .sum::<f64>();

        Some((self.t, self.y))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_euler_constant() {
        let integrator = RungeKuttaIntegrator::new_euler(0., 0., |_, _| 2., |_, y| y > &5., 1.0);
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
}
