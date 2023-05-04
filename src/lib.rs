extern crate nalgebra as na;

use na::{DMatrix, DVector, Scalar, convert, RealField, SimdValue};


/**
 * 
//     let one: T = T::one();
    // let zero: T = T::zero();
    // let max_val = T::max_value().unwrap();

    // let n = A.ncols();
    // let m = A.nrows();
    // let n_f64 : T = convert(n as f64);

    // let mut gamma_mu = DVector::<T>::from_element(n,one);
    // let mut x = DVector::<T>::from_element(n,one);
    // let mut s = DVector::<T>::from_element(n,one);
    // let mut y = DVector::<T>::zeros(m);

    // let mut x_t = DVector::<T>::from_element(n,one);
    // let mut s_t = DVector::<T>::from_element(n,one);
    // let mut y_t = DVector::<T>::zeros(m);
    
    // let mut r_primal_t = DVector::<T>::from_element(m,max_val);
    // let mut r_dual_t = DVector::<T>::from_element(n,max_val);
    // let mut r_primal = DVector::<T>::from_element(m,max_val);
    // let mut r_dual = DVector::<T>::from_element(n,max_val);
    // let mut k = 0;

    // let mut S_inv = DMatrix::<T>::zeros(n, n); // Diagonal
    // let mut X = DMatrix::<T>::zeros(n, n); // Digonal
    // let A_transpose = A.transpose();

    // let mut A_transpose_coo = CooMatrix::<T>::new(n, m);
    // let mut A_coo = CooMatrix::<T>::new(m, n);
    // for r in 0..m {
    //     for c in 0..n {
    //         let v = A[(r,c)];
    //         if v != zero {
    //             A_transpose_coo.push(c,r,v);
    //             A_coo.push(r,c,v);
    //         }

    //     }
    // }
    // let A_transpose_csc = CscMatrix::<T>::from(&A_transpose_coo);
    // let A_csc = CscMatrix::<T>::from(&A_coo);

    // let mut S_inv_coo = CooMatrix::<T>::new(n, n);
    // let mut X_coo = CooMatrix::<T>::new(n, n);

    // while k < max_it && (r_primal_t.norm() > eps || r_dual_t.norm() > eps || x.dot(&s) > eps) {
    //     //r_primal = b - A*(&x);
    //     //r_dual = c - (&A_transpose)*(&y) - (&s);

    //     r_primal_t = b - A*(&x);
    //     r_dual_t = c - (&A_transpose)*(&y) - (&s);

    //     r_primal = b - &A_csc*(&x); 
    //     r_dual = c - (&A_transpose_csc)*(&y) - (&s);
    //     let mu = x.dot(&s)/n_f64;
    //     for i in 0..n {
    //         S_inv[(i,i)] = one/s[i];
    //         S_inv_coo.push(i,i,one/s[i]);

    //         X[(i,i)] = x[i];
    //         X_coo.push(i,i,x[i]);
    //         gamma_mu[i] = gamma*mu;
    //     }

    //     let S_inv_csc = CscMatrix::<T>::from(&S_inv_coo);
    //     let X_csc = CscMatrix::<T>::from(&X_coo);

    //     let A_S_inv = A*(&S_inv);
    //     let A_S_inv_csc = (&A_csc)*(&S_inv_csc);
    //     let M_t = (&A_S_inv*(&X))*(&A_transpose); // Performance Offender
    //     let M_csc = (&A_S_inv_csc*(&X_csc))*(&A_transpose_csc); // Performance Offender
    //     let r = b + &A_S_inv*((&X)*(&r_dual_t) - (&gamma_mu)); // Performance Offender
    //     let r_csc = b + &A_S_inv_csc*((&X_csc)*(&r_dual) - (&gamma_mu));

    //     let mut M = DMatrix::<T>::zeros(M_csc.nrows(), M_csc.ncols());
    //     for (i,j,v) in M_csc.triplet_iter() {
    //         M[(i,j)] = *v;
    //     }

    //     let d_y = M_t.svd_unordered(true,true).solve(&r, svd_eps).expect("direction y failed"); // Performance Offender !!
    //     let d_s = (&r_dual_t) - (&A_transpose)*(&d_y);
    //     let d_x = -(&x) + (&S_inv)*((&gamma_mu)-(&X)*(&d_s));



    //     //let chol = CscCholesky::factor(&M_csc).expect("CSC Cholesky Failed");
    //     //let d_y_csc = chol.solve(&r_csc);
    //     let d_y_csc = M.svd_unordered(true,true).solve(&r_csc, svd_eps).expect("direction y failed");
    //     let d_s_csc = (&r_dual) - (&A_transpose_csc)*(&d_y_csc);
    //     let d_x_csc = -(&x) + (&S_inv_csc)*((&gamma_mu)-(&X_csc)*(&d_s_csc));

    //     let mut step_primal_t = max_val;
    //     let mut step_dual_t = max_val;

    //     let mut step_primal = max_val;
    //     let mut step_dual = max_val;

    //     for i in 0..n {
    //         let dx_i_t = d_x[i];
    //         let ds_i_t = d_s[i];

    //         let dx_i = d_x_csc[i];
    //         let ds_i = d_s_csc[i];

    //         if dx_i < T::zero() {
    //             let v = -x[i]/dx_i;
    //             step_primal = step_primal.min(v);
    //         }

    //         if ds_i < T::zero() {
    //             let v = -s[i]/ds_i;
    //             step_dual = step_dual.min(v);
    //         }

    //         if dx_i_t < T::zero() {
    //             let v = -x[i]/dx_i_t;
    //             step_primal_t = step_primal_t.min(v);
    //         }

    //         if ds_i_t < T::zero() {
    //             let v = -s[i]/ds_i_t;
    //             step_dual_t = step_dual_t.min(v);
    //         }
    //     }

    //     let step_max = step_primal.min(step_dual);
    //     let step = one.min(theta*step_max);

    //     x = x+d_x_csc.scale(step);
    //     y = y+d_y_csc.scale(step);
    //     s = s+d_s_csc.scale(step);

    //     let step_max_t = step_primal_t.min(step_dual_t);
    //     let step_t = one.min(theta*step_max_t);

    //     x_t = x_t+d_x.scale(step_t);
    //     y_t = y_t+d_y.scale(step_t);
    //     s_t = s_t+d_s.scale(step_t);

    //     k = k+1;
    // }
   
    // (x_t, y_t, k, r_primal_t.norm(), r_dual_t.norm())
 * 
 * 
 */

/**
 * A is Dim(MxN), x,c,s are Dim(N), y,b are Dim(M)
 * Returns (x (primal), y (dual))
 */
#[allow(non_snake_case)]
pub fn solve<T>(A: &DMatrix<T>, b: &DVector<T>, c: &DVector<T>, eps: T, theta: T, gamma: T, max_it: usize, svd_eps: T) -> (DVector<T>, DVector<T>, usize, T, T) 
    where T: Scalar + RealField + Copy + SimdValue {
    let one: T = T::one();
    let max_val = T::max_value().unwrap();
    let A_transpose = A.transpose();
    let n = A.ncols();
    let m = A.nrows();

    let mut gamma_mu = DVector::<T>::from_element(n,one);
    let mut x = DVector::<T>::from_element(n,one);
    let mut s = DVector::<T>::from_element(n,one);
    let mut y = DVector::<T>::zeros(m);
    
    let mut r_primal = DVector::<T>::from_element(m,max_val);
    let mut r_dual = DVector::<T>::from_element(n,max_val);
    let mut k = 0;

    let mut S_inv = DMatrix::<T>::zeros(n, n); // Diagonal
    let mut X = DMatrix::<T>::zeros(n, n); // Digonal
    let n_f64 : T = convert(n as f64);

    while k < max_it && (r_primal.norm() > eps || r_dual.norm() > eps || x.dot(&s) > eps) {
        r_primal = b - A*(&x);
        r_dual = c - (&A_transpose)*(&y) - (&s);
        let mu = x.dot(&s)/n_f64;
        for i in 0..n {
            S_inv[(i,i)] = one/s[i];

            X[(i,i)] = x[i];
            gamma_mu[i] = gamma*mu;
        }
        let A_S_inv= A*(&S_inv);

        let M = (&A_S_inv*(&X))*(&A_transpose); // Performance Offender
        let r = b + &A_S_inv*((&X)*(&r_dual) - (&gamma_mu)); // Performance Offender
        let d_y = M.svd_unordered(true,true).solve(&r, svd_eps).expect("direction y failed"); // Performance Offender !!
        let d_s = (&r_dual) - (&A_transpose)*(&d_y);
        let d_x = -(&x) + (&S_inv)*((&gamma_mu)-(&X)*(&d_s));

        let mut step_primal = max_val;
        let mut step_dual = max_val;

        for i in 0..n {
            let dx_i = d_x[i];
            let ds_i = d_s[i];

            if dx_i < T::zero() {
                let v = -x[i]/dx_i;
                step_primal = step_primal.min(v);
            }

            if ds_i < T::zero() {
                let v = -s[i]/ds_i;
                step_dual = step_dual.min(v);
            }
        }

        let step_max = step_primal.min(step_dual);
        let step = one.min(theta*step_max);

        x = x+d_x.scale(step);
        y = y+d_y.scale(step);
        s = s+d_s.scale(step);
        k = k+1;
    }
   
    (x, y, k, r_primal.norm(), r_dual.norm())
}