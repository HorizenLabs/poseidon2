use zkhash::{
    fields::{goldilocks::FpGoldiLocks},
    poseidon::{poseidon::Poseidon, poseidon_instance_goldilocks::{
        POSEIDON_GOLDILOCKS_8_PARAMS,
        POSEIDON_GOLDILOCKS_12_PARAMS,
        POSEIDON_GOLDILOCKS_16_PARAMS,
        POSEIDON_GOLDILOCKS_20_PARAMS,
    }},
    poseidon2::{poseidon2::Poseidon2, poseidon2_instance_goldilocks::{
        POSEIDON2_GOLDILOCKS_8_PARAMS,
        POSEIDON2_GOLDILOCKS_12_PARAMS,
        POSEIDON2_GOLDILOCKS_16_PARAMS,
        POSEIDON2_GOLDILOCKS_20_PARAMS,
    }},
    neptune::{neptune::Neptune, neptune_instances::NEPTUNE_GOLDILOCKS_PARAMS},
    gmimc::{gmimc::Gmimc, gmimc_instance_goldilocks::GMIMC_GOLDILOCKS_PARAMS},
};
type Scalar = FpGoldiLocks;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn poseidon_goldilocks(c: &mut Criterion) {
    let poseidon = Poseidon::new(&POSEIDON_GOLDILOCKS_8_PARAMS);
    let t = poseidon.get_t();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();

    c.bench_function("Poseidon Goldilocks plain (t = 12)", move |bench| {
        bench.iter(|| {
            let perm = poseidon.permutation(black_box(&input));
            black_box(perm)
        });
    });
}

fn poseidon2_goldilocks(c: &mut Criterion) {
    let poseidon2 = Poseidon2::new(&POSEIDON2_GOLDILOCKS_12_PARAMS);
    let t = poseidon2.get_t();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();

    c.bench_function("Poseidon2 Goldilocks plain (t = 12)", move |bench| {
        bench.iter(|| {
            let perm = poseidon2.permutation(black_box(&input));
            black_box(perm)
        });
    });
}

fn neptune_goldilocks(c: &mut Criterion) {
    let neptune = Neptune::new(&NEPTUNE_GOLDILOCKS_PARAMS);
    let t = neptune.get_t();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();

    c.bench_function("Neptune Goldilocks plain (t = 12)", move |bench| {
        bench.iter(|| {
            let perm = neptune.permutation(black_box(&input));
            black_box(perm)
        });
    });
}

fn gmimc_goldilocks(c: &mut Criterion) {
    let gmimc = Gmimc::new(&GMIMC_GOLDILOCKS_PARAMS);
    let t = gmimc.get_t();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();

    c.bench_function("GMiMC Goldilocks plain (t = 12)", move |bench| {
        bench.iter(|| {
            let perm = gmimc.permutation_not_opt(black_box(&input));
            black_box(perm)
        });
    });
}

fn gmimc_opt_goldilocks(c: &mut Criterion) {
    let gmimc = Gmimc::new(&GMIMC_GOLDILOCKS_PARAMS);
    let t = gmimc.get_t();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();

    c.bench_function("GMiMC (opt) Goldilocks plain (t = 12)", move |bench| {
        bench.iter(|| {
            let perm = gmimc.permutation(black_box(&input));
            black_box(perm)
        });
    });
}

fn criterion_benchmark_plain_goldilocks(c: &mut Criterion) {
    poseidon_goldilocks(c);
    poseidon2_goldilocks(c);
    neptune_goldilocks(c);
    gmimc_goldilocks(c);
    gmimc_opt_goldilocks(c);
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark_plain_goldilocks
);
criterion_main!(benches);
