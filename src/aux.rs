// use std::ops::Mul;
// use std::iter::Sum;

pub fn sum_product(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y.iter()).map(|(a, b)| a * b).sum::<f64>()
}
