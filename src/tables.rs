use ndarray::{array, Array1};

fn get_empty_array() -> Array1<f64> {
    Array1::<f64>::from_vec(vec![])
}

/// Struct storing Butcher tables or Runge-Kutta coefficients for different Runge-Kutta methods.
pub struct ButcherTableau {
    /// Number of Runge-Kutta stages
    pub num_stages: usize,

    /// Runge-Kutta order of accuracy
    pub order: usize,

    /// The `a` matrix in a Runge-Kutta Butcher table. It is two-dimensional. It is stored as a
    /// lower triangular matrix that has been flattened into a one-dimensional array slice.
    pub a: Array1<f64>,

    /// The `b` vector in a Runge-Kutta Butcher table. It is one-dimensional.
    pub b: Array1<f64>,

    /// The `c` vector in a Runge-Kutta Butcher table. It is one-dimensional.
    pub c: Array1<f64>,

    /// The `b*` coefficients in a Runge-Kutta Butcher table. These correspond
    /// to a lower order accurate scheme to be used in error-estimation and
    /// adaptive step-size control.
    pub b2: Array1<f64>, // b* in literature
}

/// Get Runge-Kutta coefficients for some method asked for.
pub fn get_rungekutta_coefficients(method_in: &str) -> Result<ButcherTableau, String> {
    let method = method_in.to_ascii_lowercase();

    match method.as_str() {
        "euler" => Ok(get_euler_coefficients()),
        "ralston" => Ok(get_ralston_coefficients()),
        "hueneuler" => Ok(get_hueneuler_coefficients()),
        "bogackishampine" => Ok(get_bogackishampine_coefficients()),
        "fehlberg" => Ok(get_fehlberg_coefficients()),
        "dormandprince" => Ok(get_dormandprince_coefficients()),
        _ => Err(format!("Runge-Kutta method -{:}- not found.", method).to_string()),
    }
}

/// Runge-Kutta coefficients table for the Euler method.
#[rustfmt::skip]
pub fn get_euler_coefficients() -> ButcherTableau {
    ButcherTableau {
        num_stages: 1,
        order: 0,
        a: get_empty_array(),
        b: array![1.],
        c: array![0.],
        b2: get_empty_array(),
    }
}

/// Runge-Kutta coefficients table for the Ralston method.
#[rustfmt::skip]
pub fn get_ralston_coefficients() -> ButcherTableau {
    ButcherTableau {
        num_stages: 2,
        order: 2,
        a: array![2./3.],
        b: array![1./4., 3./4.],
        c: array![0., 2./3.],
        b2: get_empty_array(),
    }
}

/// Runge-Kutta coefficients table for the Huen-Euler method. 
#[rustfmt::skip]
pub fn get_hueneuler_coefficients() -> ButcherTableau {
    ButcherTableau {
        num_stages: 2,
        order: 2,
        a: array![1.],
        b: array![0.5, 0.5],
        c: array![0., 1.],
        b2: array![1., 0.]
    }
}

/// Runge-Kutta coefficients table for the embedded Bogacki-Shampine method. This is a 3rd order 
/// accurate method with a 2nd order error estimator.
#[rustfmt::skip]
pub fn get_bogackishampine_coefficients() -> ButcherTableau {
    ButcherTableau {
        num_stages: 3,
        order: 3,
        a: array![
            0.5,
            0., 0.75,
            2./9., 1./3., 4./9.
            ],
        c: array![0., 0.5, 0.75, 1.],
        b: array![2./9., 1./3., 4./9., 0.],
        b2: array![7./24., 1./4., 1./3., 1./8.], 
    }
}

/// Runge-Kutta coefficients table for the embedded Fehlberg method. This is a 4th order accurate
/// method, despite the error estimator step having 5th order coefficents. 
#[rustfmt::skip]
pub fn get_fehlberg_coefficients() -> ButcherTableau {
    ButcherTableau {
        num_stages: 6,
        order: 4,
        a: array![
            0.25,
            3./32., 9./32., 
            1932./2197., -7200./2197., 7296./2197.,
            439./216., -8., 3680./513., -845./4104.,
            -8./27., 2., -3544./2565., 1859./4104., -11./40.,
            ],
        b: array![25./216., 0., 1408./2565., 2197./4104., -0.2, 0.],
        c: array![0., 0.25, 3./8., 12./13., 1., 0.5],
        b2: array![16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55.],
    }
}

/// Runge-Kutta coefficients table for the embedded Dormand-Prince method. This is a 5th order 
/// method that has an error estimator of order 4. The practice of using the higher-order method
/// to continue the integration is called "local extrapolation". 
#[rustfmt::skip]
pub fn get_dormandprince_coefficients() -> ButcherTableau {
    ButcherTableau {
        num_stages: 7,
        order: 5,
        a: array![
            0.2,
            3./40., 9./40.,
            44./45., -56./15., 32./9.,
            19372./6561., -25360./2187., 64448./6561., -212./729.,
            9017./3168., -355./33., 46732./5247., 49./176., -5103./18656.,
            35./384., 0., 500./1113., 125./192., -2187./6784., 11./84.,
            ],
        b: array![35./384., 0., 500./1113., 125./192., -2187./6784., 11./84., 0.],
        c: array![0., 0.2, 3./10., 0.8, 8./9., 1., 1.],
        b2: array![5179./57600., 0., 7571./16695., 393./640., -92097./339200., 187./2100., 1./40.],
    }
}
