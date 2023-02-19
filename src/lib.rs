extern crate nalgebra as na;

use na::{OMatrix,OVector, Scalar, Dim, DimName, default_allocator::DefaultAllocator, allocator::Allocator, convert, RealField, DimMin, DimSub, Const};


/**
 * A is Dim(MxN), x,c,s are Dim(N), y,b are Dim(M)
 * Returns (x (primal), y (dual))
 */
#[allow(non_snake_case)]
pub fn solve<T: Scalar + RealField + Copy, M: Dim + DimName + DimMin<M, Output = M> +  DimSub<Const<1>>, N: Dim + DimName + DimMin<N, Output = N>>(A: &OMatrix<T,M,N>, b: &OVector<T,M>, c: &OVector<T,N>, eps: T, max_it: usize) -> (OVector<T, N>,OVector<T,M>) 
    where  DefaultAllocator: Allocator<T, M> + Allocator<T, M, N> + Allocator<T, N, M> + Allocator<T, N, N> + Allocator<T, N> + Allocator<T, M, M> + Allocator<(T, usize), M> + Allocator<(usize, usize), M> + Allocator<T, <M as DimSub<Const<1>>>::Output>  {
    let theta : T = convert(0.95);
    let gamma: T = convert(0.9);
    
    let one: T = convert(1.0);
    let max_val = T::max_value().unwrap();
    let svd_eps: T = convert(1e-20);

    let A_transpose = A.transpose();
    let mut gamma_mu = OVector::<T,N>::from_element(one);
    let mut x = OVector::<T,N>::from_element(one);
    let mut s = OVector::<T,N>::from_element(one);
    let mut y = OVector::<T,M>::zeros();
    let mut r_primal = OVector::<T,M>::from_element(max_val);
    let mut r_dual = OVector::<T,N>::from_element(max_val);
    let mut k = 0;
    let mut S_inv = OMatrix::<T,N,N>::zeros();
    let mut X = OMatrix::<T,N,N>::zeros();
    let n : T = convert(N::USIZE as f64);

    while k < max_it && (r_primal.norm() > eps || r_dual.norm() > eps || x.dot(&s) > eps) {
        r_primal = b - A*(&x);
        r_dual = c - (&A_transpose)*(&y) - (&s);
        let mu = x.dot(&s)/n;
        for i in 0..N::USIZE {
            S_inv[(i,i)] = one/s[i];
            X[(i,i)] = x[i];
            gamma_mu[i] = gamma*mu;
        }
        let M = A*(&S_inv)*(&X)*(&A_transpose);
        let r = b + A*(&S_inv)*((&X)*(&r_dual) - (&gamma_mu));
        let d_y = M.svd(true,true).solve(&r, svd_eps).expect("direction y failed");
        let d_s = (&r_dual) - (&A_transpose)*(&d_y);
        let d_x = -(&x) + (&S_inv)*((&gamma_mu)-(&X)*(&d_s));

        let mut step_primal = max_val;
        let mut step_dual = max_val;

        for i in 0..N::USIZE {
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
   
    (x,y)
}