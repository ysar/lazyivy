use ndarray::{Array1};

/// Take the product between elements of two Rust arrays and then sum the
/// products.
pub fn sum_product(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y.iter()).map(|(a, b)| a * b).sum::<f64>()
}

/// Take the product between a float number and a one-dimensional `Array` and
/// then sum the products.
pub fn sum_product_array(x: &[f64], y: &[Array1<f64>]) -> Array1<f64> {
    let n = y[0].len();
    let mut result = Array1::<f64>::zeros(n);
    for i in 0..x.len() {
        result = result + x[i] * &y[i];
    }
    result
}