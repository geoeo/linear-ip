extern crate nalgebra as na;

use na::{OMatrix,OVector, Scalar, Dim, DimName, default_allocator::DefaultAllocator, allocator::Allocator, convert, RealField};


/**
 * A is Dim(MxN), x,c,s are Dim(N), y,b are Dim(M)
 * Returns (x (primal), y (dual))
 */
pub fn solve<T: Scalar + RealField, M: Dim + DimName, N: Dim + DimName>(A: &OMatrix<T,M,N>, b: &OVector<T,M>, c: &OVector<T,N>, eps: T, max_it: usize) -> (OVector<T, N>,OVector<T,M>) where  DefaultAllocator: Allocator<T, M> + Allocator<T, M, N>+ Allocator<T, N> {
    let theta : T = convert(0.95);
    let one: T = convert(1.0);
    let mut x = OVector::<T,N>::from_element(one.clone());
    let mut s = OVector::<T,N>::from_element(one.clone());
    let mut y = OVector::<T,M>::zeros();
    let mut r_primal = OVector::<T,M>::from_element(T::max_value().unwrap());
    let mut r_dual = OVector::<T,N>::from_element(T::max_value().unwrap());
    let mut k = 0;


   

    panic!("TODO");
}