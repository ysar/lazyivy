/// Struct storing Butcher tables or Runge-Kutta coefficients for different
/// Runge-Kutta methods that do not use a lower-order error estimator.
pub struct ButcherTableau<'t> {
    /// Number of Runge-Kutta stages
    pub s: usize,

    /// The `a` matrix in a Runge-Kutta Butcher table. It is two-dimensional.
    pub a: &'t [f64],

    /// The `b` vector in a Runge-Kutta Butcher table. It is one-dimensional.
    pub b: &'t [f64],

    /// The `c` vector in a Runge-Kutta Butcher table.
    pub c: &'t [f64],
}

/// Struct storing Butcher tables or Runge-Kutta coefficients for different
/// Runge-Kutta methods that use a lower-order error estimator. This is useful
/// for adaptive step size.
pub struct ButcherTableauAdaptive<'t> {
    /// Number of Runge-Kutta stages
    pub s: usize,

    /// Runge-Kutta order of accuracy
    pub p: usize,

    /// The `a` matrix in a Runge-Kutta Butcher table. It is two-dimensional.
    pub a: &'t [f64],

    /// The `b` vector in a Runge-Kutta Butcher table. It is one-dimensional.
    pub b: &'t [f64],

    /// The `c` vector in a Runge-Kutta Butcher table. It is one-dimensional.
    pub c: &'t [f64],

    /// The `b*` coefficients in a Runge-Kutta Butcher table. These correspond 
    /// to a lower order accurate scheme to be used in error-estimation and 
    /// adaptive step-size control.
    pub b2: &'t [f64], // b* in literature
}

/// Runge-Kutta coefficients table for the Euler method.
pub const EULER_BT: ButcherTableau = ButcherTableau {
    s: 1,
    a: &[],
    b: &[1.],
    c: &[0.],
};

/// Runge-Kutta coefficients table for the Ralston method.
#[rustfmt::skip]
pub const RALSTON_BT: ButcherTableau = ButcherTableau {
    s: 2,
    a: &[2./3.],
    b: &[1./4., 3./4.],
    c: &[0., 2./3.],
};

/// Runge-Kutta coefficients table for the adaptive Fehlberg method.
#[rustfmt::skip]
pub const FEHLBERG_BT_ADAPTIVE: ButcherTableauAdaptive = ButcherTableauAdaptive {
    s: 6,
    p: 5,
    a: &[
        0.25,
        3./32., 9./32., 
        1932./2197., -7200./2197., 7296./2197.,
        439./216., -8., 3680./513., -845./4104.,
        -8./27., 2., -3544./2565., 1859./4104., -11./40.,
        ],
    b: &[16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55.],
    c: &[0., 0.25, 3./8., 12./13., 1., 0.5],
    b2: &[25./216., 0., 1408./2565., 2197./4104., -0.2, 0.]
};

/// Runge-Kutta coefficients table for the adaptive Hueneuler method.
#[rustfmt::skip]
pub const HUENEULER_BT_ADAPTIVE: ButcherTableauAdaptive = ButcherTableauAdaptive {
    s: 2,
    p: 2,
    a: &[1.],
    b: &[0.5, 0.5],
    c: &[0., 1.],
    b2: &[1., 0.]
};
