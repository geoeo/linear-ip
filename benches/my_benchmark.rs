extern crate nalgebra as na;
use criterion::{criterion_group, criterion_main, Criterion};
use na::{Vector2,Vector1};
use linear_ip;

fn simple_f64(){
    let a = Vector2::<f64>::new(1.0,1.0);
    let b = Vector2::<f64>::new(2.0,2.0);
    let c = Vector1::<f64>::new(3.0);
    let eps: f64 = 1e-10;
    let theta : f64 = 0.95;
    let gamma : f64 = 0.1;
    let max_it = 10;

    let (_, _, _, _, _) = linear_ip::solve(&a, &b, &c, eps, theta, gamma, max_it);
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("simple_f64", |b| b.iter(|| simple_f64()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);