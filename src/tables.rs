//
// Indexing a Butcher table:
// By convention the math community uses 1-based indexing to mark the a_ij
// matrix. We use zero-based indexing, so that causes some differences and
// is a potential source for bugs.
// By our definition, the Butcher table is defined as
// ---------------------------------------------
// c0     | -
// c1     | a10, -
// c2     | a20, a21, -
// c3     | a30, a31, a32, -
// ...
// c{s-1} | a{s-1}0, a{s-1}1, ... a{s-1}{s-2}, -
// ---------------------------------------------
//        | b0, b1, b2, ..., b{s-1}
//----------------------------------------------
// For adaptive Runge-Kutta methods the lower-order coefficients are also
// provided for estimating the error
//        | b20, b21, b22, ..., b2{s-1}
//----------------------------------------------
// All coefficients are stored as 1D arrays.
// a_ij is stored as an array with indexing logic ( idx = i(i-1)/2 + j ), i.e.
// as a lower triangular matrix. This is fine for explicit Runge-Kutta methods.
//
// In the Runge-Kutta iteration, we initialize k(i=0) with the state at the
// start of the iteration. We start the stages from i=1, since a00 is undefined
// (for any explicit Runge-Kutta Method).

/// Struct storing Butcher tables or Runge-Kutta coefficients for different
/// Runge-Kutta methods that do not use a lower-order error estimator.
pub struct ButcherTableau<'t> {
    pub s: usize,
    pub a: &'t [f64],
    pub b: &'t [f64],
    pub c: &'t [f64],
}

/// Struct storing Butcher tables or Runge-Kutta coefficients for different
/// Runge-Kutta methods that use a lower-order error estimator. This is useful
/// for adaptive step size.
pub struct ButcherTableauAdaptive<'t> {
    pub s: usize,
    pub p: usize,
    pub a: &'t [f64],
    pub b: &'t [f64],
    pub c: &'t [f64],
    pub b2: &'t [f64], // b* in literature
}

pub const EULER_BT: ButcherTableau = ButcherTableau {
    s: 1,
    a: &[],
    b: &[1.],
    c: &[0.],
};

#[rustfmt::skip]
pub const RALSTON_BT: ButcherTableau = ButcherTableau {
    s: 2,
    a: &[2./3.],
    b: &[1./4., 3./4.],
    c: &[0., 2./3.],
};

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

#[rustfmt::skip]
pub const HUENEULER_BT_ADAPTIVE: ButcherTableauAdaptive = ButcherTableauAdaptive {
    s: 2,
    p: 2,
    a: &[1.],
    b: &[0.5, 0.5],
    c: &[0., 1.],
    b2: &[1., 0.]
};
