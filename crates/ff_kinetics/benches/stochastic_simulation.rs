use std::fs::File;
use std::io::{BufRead, BufReader};
use std::hint::black_box;
use criterion::criterion_group;
use criterion::criterion_main;
use criterion::Criterion;
use rand::SeedableRng;
use rand::rngs::StdRng;
use ff_structure::PairTable;
use ff_energy::ViennaRNA;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;
use ff_kinetics::LoopStructure;
use ff_kinetics::LoopStructureSSA;
use ff_kinetics::Metropolis;

fn simulate_all_from_file(path: &str) {
    let file = File::open(path).expect("Cannot open input file");
    let reader = BufReader::new(file);
    let model = ViennaRNA::default();
    let ratemodel = Metropolis::new(model.temperature(), 1.0);
    let mut rng = StdRng::seed_from_u64(42);

    let mut lines = reader.lines();
    while let Some(Ok(header)) = lines.next() {
        if !header.starts_with('>') {
            continue; // skip malformed lines
        }

        let sequence = lines.next().unwrap().unwrap();
        let structure = lines.next().unwrap().unwrap();

        let sequence = NucleotideVec::try_from(sequence.as_str()).unwrap();
        let pairings = PairTable::try_from(structure.as_str())
            .expect("invalid structure in input");
        let loops = LoopStructure::try_from((&sequence[..], &pairings, &model))
            .expect("failed to build loop structure");

        let mut simulator = LoopStructureSSA::from((loops, &ratemodel));

        simulator.simulate(
            &mut rng, 
            black_box(10.0), 
            |_t, _ti, _fl, _ls| { }
        );
    }
}

const INPUT_L50: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len50.vrna");

const INPUT_L100: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len100.vrna");

const INPUT_L250: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len250.vrna");

const INPUT_L500: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len500.vrna");

const INPUT_L750: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len750.vrna");

const INPUT_L1000: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len1000.vrna");

const INPUT_L2500: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len2500.vrna");


fn simulate_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("simulate_many");
    group.measurement_time(std::time::Duration::from_secs(50)); // increase from default 5s
    group.bench_function("simulate_l50", |b| {
        b.iter(|| simulate_all_from_file(INPUT_L50))
    });
    group.bench_function("simulate_l100", |b| {
        b.iter(|| simulate_all_from_file(INPUT_L100))
    });
    group.bench_function("simulate_l250", |b| {
        b.iter(|| simulate_all_from_file(INPUT_L250))
    });
    group.bench_function("simulate_l500", |b| {
        b.iter(|| simulate_all_from_file(INPUT_L500))
    });
    group.bench_function("simulate_l750", |b| {
        b.iter(|| simulate_all_from_file(INPUT_L750))
    });
    group.bench_function("simulate_l1000", |b| {
        b.iter(|| simulate_all_from_file(INPUT_L1000))
    });
    group.bench_function("simulate_l2500", |b| {
        b.iter(|| simulate_all_from_file(INPUT_L2500))
    });
    group.finish();
}

criterion_group!(benches, simulate_benchmark);
criterion_main!(benches);

