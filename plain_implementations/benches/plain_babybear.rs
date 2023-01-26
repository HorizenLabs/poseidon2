use zkhash::{
    fields::{babybear::FpBabyBear},
    poseidon::{poseidon::Poseidon, poseidon_instance_babybear::POSEIDON_BABYBEAR_PARAMS},
    poseidon2::{poseidon2::Poseidon2, poseidon2_instance_babybear::POSEIDON2_BABYBEAR_PARAMS},
    neptune::{neptune::Neptune, neptune_instances::NEPTUNE_BABYBEAR_PARAMS},
    gmimc::{gmimc::Gmimc, gmimc_instance_babybear::GMIMC_BABYBEAR_PARAMS},
};
type Scalar = FpBabyBear;
// use ark_ff::UniformRand;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn poseidon_babybear(c: &mut Criterion) {
    let poseidon = Poseidon::new(&POSEIDON_BABYBEAR_PARAMS);
    let t = poseidon.get_t();
    // let input: Vec<Scalar> = (0..t).map(|_| Scalar::rand(&mut ark_std::rand::thread_rng())).collect();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();
    c.bench_function("Poseidon BabyBear plain (t = 24)", move |bench| {
        bench.iter(|| {
            let perm = poseidon.permutation(black_box(&input));
            black_box(perm)
        });
    });
}

fn poseidon2_babybear(c: &mut Criterion) {
    let poseidon2 = Poseidon2::new(&POSEIDON2_BABYBEAR_PARAMS);
    let t = poseidon2.get_t();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();

    c.bench_function("Poseidon2 BabyBear plain (t = 24)", move |bench| {
        bench.iter(|| {
            let perm = poseidon2.permutation(black_box(&input));
            black_box(perm)
        });
    });
}

fn neptune_babybear(c: &mut Criterion) {
    let neptune = Neptune::new(&NEPTUNE_BABYBEAR_PARAMS);
    let t = neptune.get_t();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();

    c.bench_function("Neptune BabyBear plain (t = 24)", move |bench| {
        bench.iter(|| {
            let perm = neptune.permutation(black_box(&input));
            black_box(perm)
        });
    });
}

fn gmimc_babybear(c: &mut Criterion) {
    let gmimc = Gmimc::new(&GMIMC_BABYBEAR_PARAMS);
    let t = gmimc.get_t();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();

    c.bench_function("GMiMC BabyBear plain (t = 24)", move |bench| {
        bench.iter(|| {
            let perm = gmimc.permutation_not_opt(black_box(&input));
            black_box(perm)
        });
    });
}

fn gmimc_opt_babybear(c: &mut Criterion) {
    let gmimc = Gmimc::new(&GMIMC_BABYBEAR_PARAMS);
    let t = gmimc.get_t();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();

    c.bench_function("GMiMC (opt) BabyBear plain (t = 24)", move |bench| {
        bench.iter(|| {
            let perm = gmimc.permutation(black_box(&input));
            black_box(perm)
        });
    });
}

fn criterion_benchmark_plain_babybear(c: &mut Criterion) {
    poseidon_babybear(c);
    poseidon2_babybear(c);
    neptune_babybear(c);
    gmimc_babybear(c);
    gmimc_opt_babybear(c);
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark_plain_babybear
);
criterion_main!(benches);
