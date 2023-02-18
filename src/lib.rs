extern crate nalgebra as na;

use na::{OMatrix,OVector, Scalar, Dim, DimName, default_allocator::DefaultAllocator, allocator::Allocator, convert, RealField, DimMin};


/**
 * A is Dim(MxN), x,c,s are Dim(N), y,b are Dim(M)
 * Returns (x (primal), y (dual))
 */
#[allow(non_snake_case)]
pub fn solve<T: Scalar + RealField + Copy, M: Dim + DimName, N: Dim + DimName + DimMin<N, Output = N>>(A: &OMatrix<T,M,N>, b: &OVector<T,M>, c: &OVector<T,N>, eps: T, max_it: usize) -> (OVector<T, N>,OVector<T,M>) 
    where  DefaultAllocator: Allocator<T, M> + Allocator<T, M, N> + Allocator<T, N, M> + Allocator<T, N, N> + Allocator<T, N> + Allocator<T, M, M> {
    let theta : T = convert(0.95);
    let one: T = convert(1.0);
    let A_tranpose = A.transpose();
    let mut x = OVector::<T,N>::from_element(one);
    let mut s = OVector::<T,N>::from_element(one);
    let mut y = OVector::<T,M>::zeros();
    let mut r_primal = OVector::<T,M>::from_element(T::max_value().unwrap());
    let mut r_dual = OVector::<T,N>::from_element(T::max_value().unwrap());
    let mut k = 0;
    let mut S_inv = OMatrix::<T,N,N>::zeros();
    let mut X = OMatrix::<T,N,N>::zeros();
    let n : T = convert(N::USIZE as f64);

    while k < max_it && (r_primal.norm() > eps || r_dual.norm() > eps || x.dot(&s) > eps) {
        r_primal = b - A*(&x);
        r_dual = c - (&A_tranpose)*(&y) - (&s);
        let mu = x.dot(&s)/n;
        for i in 0..N::USIZE {
            S_inv[(i,i)] = one/s[i];
            X[(i,i)] = x[i];
        }
        let M = A*(&S_inv)*(&X)*(&A_tranpose);

       

    }
   

    panic!("TODO");
}