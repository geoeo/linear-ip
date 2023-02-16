extern crate nalgebra as na;

use na::{OMatrix,OVector, Scalar, Dim, default_allocator::DefaultAllocator, allocator::Allocator};


/**
 * A is Dim(MxN), x,c are Dim(N), y,b are Dim(M)
 * Returns (x (primal), y (dual))
 */
pub fn solve<T: Scalar, M: Dim, N: Dim>(A: &OMatrix<T,M,N>, b: &OVector<T,M>, c: &OVector<T,N>) -> (OVector<T, N>,OVector<T,M>) where  DefaultAllocator: Allocator<T, M> + Allocator<T, M, N>+ Allocator<T, N> {
    panic!("TODO");
}