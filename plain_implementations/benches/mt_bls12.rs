use criterion::{black_box, criterion_group, criterion_main, Criterion};
use zkhash::{
    fields::{bls12::FpBLS12},
    merkle_tree::merkle_tree_fp::MerkleTree,
    neptune::{neptune::Neptune, neptune_instances::NEPTUNE_BLS_PARAMS},
    poseidon::{poseidon::Poseidon, poseidon_instance_bls12::POSEIDON_BLS_PARAMS},
};
type Scalar = FpBLS12;

fn sample_set(set_size: usize) -> Vec<Scalar> {
    (0..set_size).map(|i| Scalar::from(i as u64)).collect()
}

fn poseidon(c: &mut Criterion, log_set_size: usize) {
    let perm = Poseidon::new(&POSEIDON_BLS_PARAMS);
    let mut mt = MerkleTree::new(perm);
    let set_size = 1 << log_set_size;
    let set: Vec<Scalar> = sample_set(set_size);

    let id = format!("Poseidon BLS12 MT (set_size = 2^{})", log_set_size);

    c.bench_function(&id, move |bench| {
        bench.iter(|| {
            mt.accumulate(black_box(&set));
        });
    });
}

fn neptune(c: &mut Criterion, log_set_size: usize) {
    let perm = Neptune::new(&NEPTUNE_BLS_PARAMS);
    let mut mt = MerkleTree::new(perm);
    let set_size = 1 << log_set_size;
    let set: Vec<Scalar> = sample_set(set_size);

    let id = format!("Neptune BLS12 MT (set_size = 2^{})", log_set_size);

    c.bench_function(&id, move |bench| {
        bench.iter(|| {
            mt.accumulate(black_box(&set));
        });
    });
}

fn criterion_benchmark_mt_bls(c: &mut Criterion) {
    let log_set_sizes = vec![20];

    for log_set_size in log_set_sizes {
        poseidon(c, log_set_size);
        neptune(c, log_set_size);
    }
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = criterion_benchmark_mt_bls
);
criterion_main!(benches);
