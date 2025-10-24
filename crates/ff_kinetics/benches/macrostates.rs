use criterion::Criterion;
use criterion::criterion_group;
use criterion::criterion_main;

use ff_structure::DotBracketVec;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;
use ff_energy::ViennaRNA;

use ff_kinetics::macrostates::Macrostate;
use ff_kinetics::Metropolis;

pub fn macrostate_neighbors(c: &mut Criterion) {
    let mut group = c.benchmark_group("Macrostates");
    //group.measurement_time(std::time::Duration::from_secs(50)); // increase from default 5s
    let seq = NucleotideVec::from_lossy("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC");
    let energy_model = ViennaRNA::default();
    let rate_model = Metropolis::new(energy_model.temperature(), 1.0);

    let s1 = DotBracketVec::try_from(".((((....)))).((((........))))...............").unwrap();
    let s2 = DotBracketVec::try_from(".((((....)))).((((.(....).))))...............").unwrap();
    let s3 = DotBracketVec::try_from(".((((....))))..(((........)))................").unwrap();
    let s4 = DotBracketVec::try_from(".((((....)))).((((.(.....)))))...............").unwrap();
    let s5 = DotBracketVec::try_from(".(((......))).((((........))))...............").unwrap();
    let s6 = DotBracketVec::try_from("..(((....)))..((((........))))...............").unwrap();
    let s7 = DotBracketVec::try_from(".(((......)))..(((........)))................").unwrap();
    let s8 = DotBracketVec::try_from(".(((.(...)))).((((........))))...............").unwrap();

    let macrostate = Macrostate::from_list(
        "bench_state",
        &seq,
        &[s1, s2, s3, s4, s5, s6, s7, s8], 
        &energy_model,
    );

    group.bench_function("Neighborhood generation", |b| {
        b.iter(|| {
            let _ = Macrostate::from_macrostate_neighbors(
                &macrostate, 
                &seq, 
                &energy_model, 
                &rate_model
            );
        });
    });
}

criterion_group!(benches, macrostate_neighbors);
criterion_main!(benches);

