extern crate nalgebra as na;
extern crate nalgebra_sparse;

use na::{DMatrix, DVector, Scalar, convert, RealField, SimdValue};
use nalgebra_sparse::{CooMatrix,CscMatrix};

/**
 * A is Dim(MxN), x,c,s are Dim(N), y,b are Dim(M)
 * Returns (x (primal), y (dual))
 */
#[allow(non_snake_case)]
pub fn solve<T>(A: &DMatrix<T>, b: &DVector<T>, c: &DVector<T>, eps: T, theta: T, gamma: T, max_it: usize, svd_eps: T) -> (DVector<T>, DVector<T>, usize, T, T) 
    where T: Scalar + RealField + Copy + SimdValue {
    let one = T::one();
    let zero = T::zero();
    let max_val = T::max_value().unwrap();
    let n = A.ncols();
    let m = A.nrows();

    let mut gamma_mu = DVector::<T>::from_element(n,one);

    let mut x_sp = DVector::<T>::from_element(n,one);
    let mut s_sp = DVector::<T>::from_element(n,one);
    let mut y_sp = DVector::<T>::zeros(m);

    let mut r_primal_sp = DVector::<T>::from_element(m,max_val);
    let mut r_dual_sp = DVector::<T>::from_element(n,max_val);

    let mut k = 0;

    let mut A_transpose_coo = CooMatrix::<T>::new(n, m);
    let mut A_coo = CooMatrix::<T>::new(m, n);
    for r in 0..m {
        for c in 0..n {
            let v = A[(r,c)];
            if v != zero {
                A_transpose_coo.push(c,r,v);
                A_coo.push(r,c,v);
            }

        }
    }
    let A_transpose_sp = CscMatrix::<T>::from(&A_transpose_coo);
    let A_sp = CscMatrix::<T>::from(&A_coo);

    let mut S_inv_coo = CooMatrix::<T>::new(n, n);
    let mut X_coo = CooMatrix::<T>::new(n, n);

    let n_f64 : T = convert(n as f64);

    while k < max_it && (r_primal_sp.norm() > eps || r_dual_sp.norm() > eps || x_sp.dot(&s_sp) > eps) {

        r_primal_sp = b - &A_sp*(&x_sp);
        r_dual_sp = c - (&A_transpose_sp)*(&y_sp) - (&s_sp);

        let mu = x_sp.dot(&s_sp)/n_f64;

        for i in 0..n {
            let v_s = one/s_sp[i];
            if v_s != zero {
                S_inv_coo.push(i,i,one/s_sp[i]);
            }

            let v_x = x_sp[i];
            if v_x != zero {
                X_coo.push(i,i,v_x);
            }
            gamma_mu[i] = gamma*mu;
        }

        let S_inv_sp = CscMatrix::<T>::from(&S_inv_coo);
        let X_sp = CscMatrix::<T>::from(&X_coo);
        let A_S_inv_sp = (&A_sp)*(&S_inv_sp);

        let M_sp = (&A_S_inv_sp*(&X_sp))*(&A_transpose_sp); // Performance Offender
        let r_sp = b + &A_S_inv_sp*((&X_sp)*(&r_dual_sp) - (&gamma_mu));

        let mut M_sp_d = DMatrix::<T>::zeros(m, m);
        for (i,j,v) in M_sp.triplet_iter() {
            M_sp_d[(i,j)] = *v;
        }

        let d_y_sp = M_sp_d.svd_unordered(true,true).solve(&r_sp, svd_eps).expect("direction y failed");
        let d_s_sp = (&r_dual_sp) - (&A_transpose_sp)*(&d_y_sp);
        let d_x_sp = -(&x_sp) + (&S_inv_sp)*((&gamma_mu)-(&X_sp)*(&d_s_sp));

        let mut step_primal_sp = max_val;
        let mut step_dual_sp = max_val;

        for i in 0..n {

            let dx_i_sp = d_x_sp[i];
            let ds_i_sp = d_s_sp[i];

            if dx_i_sp < T::zero() {
                let v = -x_sp[i]/dx_i_sp;
                  step_primal_sp = step_primal_sp.min(v);
            }

            if ds_i_sp < T::zero() {
                let v = -s_sp[i]/ds_i_sp;
                step_dual_sp = step_dual_sp.min(v);
            }
        }

        let step_max_sp = step_primal_sp.min(step_dual_sp);
        let step_sp = one.min(theta*step_max_sp);

        x_sp = x_sp+d_x_sp.scale(step_sp);
        y_sp = y_sp+d_y_sp.scale(step_sp);
        s_sp = s_sp+d_s_sp.scale(step_sp);

        S_inv_coo.clear_triplets();
        X_coo.clear_triplets();

        k = k+1;
    }
   
    (x_sp, y_sp, k, r_primal_sp.norm(), r_dual_sp.norm())
}