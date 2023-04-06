extern crate nalgebra as na;
#[macro_use]
extern crate approx;

use na::{Vector2,Vector1,DVector, DMatrix};
use linear_ip;

/**
 * (Dual) Max 2*x_1 + 2*x_2 s.t. x_1 + x_2 <= 3 and x_1, x_2 >= 0
 */
#[test]
fn simple_f64(){
    let a = Vector2::<f64>::new(1.0,1.0);
    let b = Vector2::<f64>::new(2.0,2.0);
    let c = Vector1::<f64>::new(3.0);
    let eps: f64 = 1e-10;
    let theta : f64 = 0.95;
    let gamma : f64 = 0.1;
    let max_it = 10;

    let (_, y, _, _, _) = linear_ip::solve(&a, &b, &c, eps, theta, gamma, max_it);
    let y1 = y[0];
    let y2 = y[1];
    assert_relative_eq!(y1,1.5f64, epsilon  = eps);
    assert_relative_eq!(y2,1.5f64, epsilon  = eps);
}

/**
 * (Dual) Max 2*x_1 + 2*x_2 s.t. x_1 + x_2 <= 3 and x_1, x_2 >= 0
 */
#[test]
fn simple_f64_dyn(){
    let a = DMatrix::<f64>::from_vec(2,1, vec![1.0,1.0]);
    let b = DVector::<f64>::from_vec(vec![2.0,2.0]);
    let c = DVector::<f64>::from_vec(vec![3.0]);
    let eps: f64 = 1e-10;
    let theta : f64 = 0.95;
    let gamma : f64 = 0.1;
    let max_it = 10;

    let (_, y, _, _, _) = linear_ip::solve(&a, &b, &c, eps, theta, gamma, max_it);
    let y1 = y[0];
    let y2 = y[1];
    assert_relative_eq!(y1,1.5f64, epsilon  = eps);
    assert_relative_eq!(y2,1.5f64, epsilon  = eps);
}

/**
 * (Dual) Max 2*x_1 + 2*x_2 s.t. x_1 + x_2 <= 3 and x_1, x_2 >= 0
 */
#[test]
fn simple_f32(){
    let a = Vector2::<f32>::new(1.0,1.0);
    let b = Vector2::<f32>::new(2.0,2.0);
    let c = Vector1::<f32>::new(3.0);
    let eps: f32 = 1e-10;
    let theta : f32 = 0.95;
    let gamma : f32 = 0.1;
    let max_it = 10;

    let (_, y, _, _, _) = linear_ip::solve(&a, &b, &c, eps, theta, gamma, max_it);
    let y1 = y[0];
    let y2 = y[1];
    assert_relative_eq!(y1,1.5f32, epsilon  = eps);
    assert_relative_eq!(y2,1.5f32, epsilon  = eps);
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
    let theta : f32 = 0.95;
    let gamma : f32 = 0.1;
    let max_it = 10;

    let (_, y, _, _, _) = linear_ip::solve(&a, &b, &c, eps, theta, gamma, max_it);
    let y1 = y[0];
    let y2 = y[1];
    assert_relative_eq!(y1,1.5f32, epsilon  = eps);
    assert_relative_eq!(y2,1.5f32, epsilon  = eps);
}
