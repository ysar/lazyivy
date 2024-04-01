// Linalg routines
use lazyivy::linalg::dot_product;
use lazyivy::linalg::sum_product;
use lazyivy::linalg::sum_product_array;
use lazyivy::linalg::vector_matrix_product;

//
use ndarray::{Array1, Array2};
use ndarray_rand::rand_distr::Uniform;
use ndarray_rand::RandomExt;

//
use std::time::{Instant};

// Array generation
fn gen_1d_array(size: usize) -> Array1<f64> {
    let random_array =
      Array1::random(size, Uniform::new(0.0, 1.0));
    random_array
}

fn gen_2d_array(nb_rows: usize, nb_columns: usize) -> Array2<f64> {
    let random_array =
      Array2::random((nb_rows, nb_columns), Uniform::new(0.0, 1.0));
    
    random_array
}

fn main() {

    //
    const ARRAY_SIZE: usize = 10_000;

    // Generate two random arrays

    //  -> 1D vectors
    let vector1: Array1<f64> = gen_1d_array(ARRAY_SIZE);
    let vector2: Array1<f64> = gen_1d_array(ARRAY_SIZE);

    //  -> 2D matrix
    let matrix: Array2<f64> =
      gen_2d_array(ARRAY_SIZE, ARRAY_SIZE);

    // Dot product timings
    let timer = Instant::now();
    let out1 = sum_product(&vector1, &vector2);
    let time1 = timer.elapsed();

    let timer = Instant::now();
    let out2 = dot_product(&vector1, &vector2);
    let time2 = timer.elapsed();

    println!("Sum product: {:?} in {:?}", out1, time1);
    println!("Dot product: {:?} in {:?}", out2, time2);

    println!("");

    // Row-Matrix product
    let timer = Instant::now();
    let _product1: Array1<f64> = sum_product_array(&vector1, &matrix);
    let time3 = timer.elapsed();

    let timer = Instant::now();
    let _product2: Array1<f64> = vector_matrix_product(&vector1, &matrix);
    let time4 = timer.elapsed();

    println!("Sum product array: {:?}", time3);
    println!("Vector matrix product: {:?}", time4);
    
}