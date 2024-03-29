extern crate nalgebra as na;
#[macro_use]
extern crate approx;

use na::{DVector, DMatrix};
use linear_ip;


/**
 * (Dual) Max 2*x_1 + 2*x_2 s.t. x_1 + x_2 <= 3 and x_1, x_2 >= 0
 */
#[test]
fn simple_f64_dyn(){
    let a = DMatrix::<f64>::from_vec(2,1, vec![1.0,1.0]);
    let b = DVector::<f64>::from_vec(vec![2.0,2.0]);
    let c = DVector::<f64>::from_vec(vec![3.0]);
    let eps: f64 = 1e-10;
    let theta : f64 = 0.25;
    let gamma : f64 = 0.5;
    let max_it = 1000;

    let (_, y, _, _, _) = linear_ip::solve(&a, &b, &c, eps, theta, gamma, max_it, 1e-20);
    let y1 = y[0];
    let y2 = y[1];
    assert_relative_eq!(y1,1.5f64, epsilon  = eps);
    assert_relative_eq!(y2,1.5f64, epsilon  = eps);
}

/**
 * (Dual) Max 2*x_1 + 2*x_2 s.t. x_1 + x_2 <= 3 and x_1, x_2 >= 0
 */
#[test]
fn simple_f32_dyn(){
    let a = DMatrix::<f32>::from_vec(2,1,vec![1.0,1.0]);
    let b = DVector::<f32>::from_vec(vec![2.0,2.0]);
    let c = DVector::<f32>::from_vec(vec![3.0]);
    let eps: f32 = 1e-10;
    let theta : f32 = 0.25;
    let gamma : f32 = 0.5;
    let max_it = 100;

    let (_, y, _, _, _) = linear_ip::solve(&a, &b, &c, eps, theta, gamma, max_it, 1e-20);
    let y1 = y[0];
    let y2 = y[1];
    assert_relative_eq!(y1,1.5f32, epsilon  = eps);
    assert_relative_eq!(y2,1.5f32, epsilon  = eps);
}
