/// Struct for handling errors that arise while building the `RungeKutta` struct.
#[derive(Debug)]
pub enum BuilderError {
    /// Error type for mismatch between array sizes.
    InconsistentSize(String),

    /// Error type for methods for which adaptive step-size has not been implemented.
    AdaptiveNotImplemented,
}

impl std::fmt::Display for BuilderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InconsistentSize(s) => write!(
                f,
                "Array size of {s} does not match size of initial condition."
            ),
            Self::AdaptiveNotImplemented => write!(
                f,
                "Adaptive step-size not implemented for selected Runge-Kutta method."
            ),
        }
    }
}
