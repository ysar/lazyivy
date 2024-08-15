use thiserror::Error;

/// Struct for handling errors that arise while building the `RungeKutta` struct.
#[derive(Error, Debug)]
pub enum BuilderError {
    
    /// Error type for mismatch between array sizes.
    #[error("Array size of `{0}` does not match size of initial condition.")]
    InconsistentSize(String),
    
    /// Error type for methods for which adaptive step-size has not been implemented.
    #[error("Adaptive step-size not implemented for method : {0}")]
    AdaptiveNotImplemented(String),
    
    /// Error type for an unknown Runge-Kutta method.
    #[error("Unknown Runge-Kutta method: {0}")]
    UnknownRungeKuttaMethod(String),
}
