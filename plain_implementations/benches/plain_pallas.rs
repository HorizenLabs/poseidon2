use zkhash::{
    fields::{pallas::FpPallas},
    poseidon::{poseidon::Poseidon, poseidon_instance_pallas::{
        POSEIDON_PALLAS_3_PARAMS,
        POSEIDON_PALLAS_4_PARAMS,
        POSEIDON_PALLAS_8_PARAMS
    }},
    poseidon2::{poseidon2::Poseidon2, poseidon2_instance_pallas::POSEIDON2_PALLAS_PARAMS},
    gmimc::{gmimc::Gmimc, gmimc_instance_pallas::GMIMC_PALLAS_3_PARAMS},
};
type Scalar = FpPallas;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn poseidon_pallas(c: &mut Criterion) {
    let poseidon = Poseidon::new(&POSEIDON_PALLAS_3_PARAMS);
    let t = poseidon.get_t();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();

    c.bench_function("Poseidon Pallas plain (t = 3)", move |bench| {
        bench.iter(|| {
            let perm = poseidon.permutation(black_box(&input));
            black_box(perm)
        });
    });
}

fn poseidon2_pallas(c: &mut Criterion) {
    let poseidon2 = Poseidon2::new(&POSEIDON2_PALLAS_PARAMS);
    let t = poseidon2.get_t();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();

    c.bench_function("Poseidon2 Pallas plain (t = 3)", move |bench| {
        bench.iter(|| {
            let perm = poseidon2.permutation(black_box(&input));
            black_box(perm)
        });
    });
}

fn gmimc_pallas(c: &mut Criterion) {
    let gmimc = Gmimc::new(&GMIMC_PALLAS_3_PARAMS);
    let t = gmimc.get_t();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();

    c.bench_function("GMiMC Pallas plain (t = 3)", move |bench| {
        bench.iter(|| {
            let perm = gmimc.permutation_not_opt(black_box(&input));
            black_box(perm)
        });
    });
}

fn gmimc_opt_pallas(c: &mut Criterion) {
    let gmimc = Gmimc::new(&GMIMC_PALLAS_3_PARAMS);
    let t = gmimc.get_t();
    let input: Vec<Scalar> = (0..t).map(|i| Scalar::from(i as u64)).collect();

    c.bench_function("GMiMC (opt) Pallas plain (t = 3)", move |bench| {
        bench.iter(|| {
            let perm = gmimc.permutation(black_box(&input));
            black_box(perm)
        });
    });
}

fn criterion_benchmark_plain_pallas(c: &mut Criterion) {
    poseidon_pallas(c);
    poseidon2_pallas(c);
    gmimc_pallas(c);
    gmimc_opt_pallas(c);
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark_plain_pallas
);
criterion_main!(benches);
