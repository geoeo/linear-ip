extern crate nalgebra as na;
use criterion::{criterion_group, criterion_main, Criterion};
use na::{DMatrix,DVector};
use linear_ip;

fn simple_f64(){
    let a = DMatrix::<f64>::from_vec(2,1,vec![1.0,1.0]);
    let b = DVector::<f64>::from_vec(vec![2.0,2.0]);
    let c = DVector::<f64>::from_vec(vec![3.0]);

    let eps: f64 = 1e-10;
    let theta : f64 = 0.95;
    let gamma : f64 = 0.1;
    let max_it = 10;

    let (_, _, _, _, _) = linear_ip::solve(&a, &b, &c, eps, theta, gamma, max_it, 1e-20);
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("simple_f64", |b| b.iter(|| simple_f64()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);