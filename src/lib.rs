extern crate nalgebra as na;

use na::{OMatrix, OVector, Vector, Matrix, Scalar, Dim, default_allocator::DefaultAllocator, allocator::Allocator, convert, RealField, DimMin, DimSub, Const, storage::Owned, SimdValue};


/**
 * A is Dim(MxN), x,c,s are Dim(N), y,b are Dim(M)
 * Returns (x (primal), y (dual))
 */
#[allow(non_snake_case)]
pub fn solve<'a, T, M, N>(A: &'a OMatrix<T,M,N>, b: &'a OVector<T,M>, c: &'a OVector<T,N>, eps: T, theta: T, gamma: T, max_it: usize, svd_eps: T) -> (OVector<T, N>,OVector<T,M>, usize, T, T) 
    where 
        T: Scalar + RealField + Copy + SimdValue, 
        M: Dim + DimMin<M, Output = M> + DimSub<Const<1>>, 
        N: Dim + DimMin<N, Output = N>,
        DefaultAllocator: Allocator<T, Const<1>, N> + Allocator<T, M> + Allocator<T, M, N> + Allocator<T, N, M> + Allocator<T, N, N> + Allocator<T, N> + Allocator<T, M, M> + Allocator<(T, usize), M> 
        + Allocator<(usize, usize), M> + Allocator<T, <M as DimSub<Const<1>>>::Output> {
    let one: T = convert(1.0);
    let max_val = T::max_value().unwrap();
    let A_transpose = A.transpose();
    let n = A.ncols();
    let m = A.nrows();
    let dim_n = Dim::from_usize(n);
    let dim_m = Dim::from_usize(m);

    let mut gamma_mu = Vector::<T,N, Owned<T,N,Const<1>>>::from_element_generic(dim_n, Const::<1>,one);
    let mut x = Vector::<T,N, Owned<T,N,Const<1>>>::from_element_generic(dim_n, Const::<1>,one);
    let mut s = Vector::<T,N, Owned<T,N,Const<1>>>::from_element_generic(dim_n, Const::<1>,one);
    let mut y = Vector::<T,M, Owned<T,M,Const<1>>>::zeros_generic(dim_m, Const::<1>);
    let mut r_primal =  Vector::<T,M, Owned<T,M,Const<1>>>::from_element_generic(Dim::from_usize(m), Const::<1>,max_val);
    let mut r_dual = Vector::<T,N, Owned<T,N,Const<1>>>::from_element_generic(dim_n, Const::<1>,max_val);
    let mut k = 0;
    let mut S_inv = Matrix::<T,N,N,Owned<T,N,N>>::zeros_generic(dim_n, dim_n); // Diagonal
    let mut X = Matrix::<T,N,N,Owned<T,N,N>>::zeros_generic(dim_n, dim_n); // Digonal
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