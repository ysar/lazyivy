/// Struct storing Butcher tables or Runge-Kutta coefficients for different Runge-Kutta methods.
pub struct ButcherTableau<'t> {
    /// Number of Runge-Kutta stages
    pub s: usize,

    /// Runge-Kutta order of accuracy
    pub p: usize,

    /// The `a` matrix in a Runge-Kutta Butcher table. It is two-dimensional. It is stored as a
    /// lower triangular matrix that has been flattened into a one-dimensional array slice.
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
    p: 0,
    a: &[],
    b: &[1.],
    c: &[0.],
    b2: &[],
};

/// Runge-Kutta coefficients table for the Ralston method.
#[rustfmt::skip]
pub const RALSTON_BT: ButcherTableau = ButcherTableau {
    s: 2,
    p: 2,
    a: &[2./3.],
    b: &[1./4., 3./4.],
    c: &[0., 2./3.],
    b2: &[],
};

/// Runge-Kutta coefficients table for the Huen-Euler method. This is a second order 
/// method with an error estimate of order 1.
#[rustfmt::skip]
pub const HUENEULER_BT: ButcherTableau = ButcherTableau {
    s: 2,
    p: 2,
    a: &[1.],
    b: &[0.5, 0.5],
    c: &[0., 1.],
    b2: &[1., 0.]
};

/// Runge-Kutta coefficients table for the embedded Bogacki-Shampine method. This is a 3rd order 
/// accurate method with a 2nd order error estimator.
#[rustfmt::skip]
pub const BOGACKISHAMPINE_BT: ButcherTableau = ButcherTableau {
    s: 3,
    p: 3,
    a: &[
        0.5,
        0., 0.75,
        2./9., 1./3., 4./9.
        ],
    c: &[0., 0.5, 0.75, 1.],
    b: &[2./9., 1./3., 4./9., 0.],
    b2: &[7./24., 1./4., 1./3., 1./8.], 
};

/// Runge-Kutta coefficients table for the embedded Fehlberg method. This is a 4th order accurate
/// method, despite the error estimator step having 5th order coefficents. 
#[rustfmt::skip]
pub const FEHLBERG_BT: ButcherTableau = ButcherTableau {
    s: 6,
    p: 4,
    a: &[
        0.25,
        3./32., 9./32., 
        1932./2197., -7200./2197., 7296./2197.,
        439./216., -8., 3680./513., -845./4104.,
        -8./27., 2., -3544./2565., 1859./4104., -11./40.,
        ],
    b: &[25./216., 0., 1408./2565., 2197./4104., -0.2, 0.],
    c: &[0., 0.25, 3./8., 12./13., 1., 0.5],
    b2: &[16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55.],
};

/// Runge-Kutta coefficients table for the embedded Dormand-Prince method. This is a 5th order 
/// method that has an error estimator of order 4. The practice of using the higher-order method
/// to continue the integration is called "local extrapolation". 
#[rustfmt::skip]
pub const DORMANDPRINCE_BT: ButcherTableau = ButcherTableau {
    s: 7,
    p: 5,
    a: &[
        0.2,
        3./40., 9./40.,
        44./45., -56./15., 32./9.,
        19372./6561., -25360./2187., 64448./6561., -212./729.,
        9017./3168., -355./33., 46732./5247., 49./176., -5103./18656.,
        35./384., 0., 500./1113., 125./192., -2187./6784., 11./84.,
        ],
    b: &[35./384., 0., 500./1113., 125./192., -2187./6784., 11./84., 0.],
    c: &[0., 0.2, 3./10., 0.8, 8./9., 1., 1.],
    b2: &[5179./57600., 0., 7571./16695., 393./640., -92097./339200., 187./2100., 1./40.],
};
