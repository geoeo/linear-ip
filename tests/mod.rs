extern crate nalgebra as na;
#[macro_use]
extern crate approx;

use na::{Vector2,Vector1};
use linear_ip;

/**
 * (Dual) Max 2*x_1 + 2*x_2 s.t. x_1 + x_2 <= 3 and x_1, x_2 >= 0
 */
#[test]
fn simple(){
    let a = Vector2::<f64>::new(1.0,1.0);
    let b = Vector2::<f64>::new(2.0,2.0);
    let c = Vector1::<f64>::new(3.0);
    let eps: f64 = 1e-5;
    let max_it = 100;

    let (_,y) = linear_ip::solve(&a, &b, &c, eps, max_it);
    let y1 = y[0];
    let y2 = y[1];
    assert_relative_eq!(y1,1.5f64, epsilon  = eps);
    assert_relative_eq!(y2,1.5f64, epsilon  = eps);
}
