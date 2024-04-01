use ndarray::Array1;
use ndarray::Array2;

/// Take the product between elements of two Rust arrays and then sum the products.
pub fn sum_product(x: &Array1<f64>, y: &Array1<f64>) -> f64 {
    x.iter().zip(y.iter()).map(|(a, b)| a * b).sum::<f64>()
}

/// Dot product with an explicit for loop
pub fn dot_product(x: &Array1<f64>, y: &Array1<f64>) -> f64 {
    x.dot(y)
}

/// Take the product between a float number and a one-dimensional `Array` and then sum the products.
pub fn sum_product_array(x: &Array1<f64>, y: &Array2<f64>) -> Array1<f64> {
    let (nb_rows, nb_cols) = y.dim();

    let mut result =
      Array1::<f64>::zeros(nb_cols);
    for row_id in 0..nb_rows {
        result = result + x[row_id] * &y.row(row_id);
    }
    result
}

/// Multiply a matrix on the left with a vector
pub fn vector_matrix_product(vector: &Array1<f64>, matrix: &Array2<f64>) -> Array1<f64> {
    let (nb_rows, nb_cols) = matrix.dim();

    let mut result = 
      Array1::<f64>::zeros(nb_cols);

    for row_id in 0..nb_rows {
        let factor: f64 = vector[row_id];

        for col_id in 0..nb_cols {
           result[col_id] += factor * matrix[[row_id, col_id]];
        }
    }

    result
}