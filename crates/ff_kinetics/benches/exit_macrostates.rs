use criterion::Criterion;
use criterion::criterion_group;
use criterion::criterion_main;

use ff_structure::DotBracketVec;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;
use ff_energy::ViennaRNA;

use ff_kinetics::Metropolis;
use ff_kinetics::Macrostate;
use ff_kinetics::commit_and_delay::ExitMacrostate;

pub fn exit_macrostates(c: &mut Criterion) {
    let mut group = c.benchmark_group("Macrostates");
    //group.measurement_time(std::time::Duration::from_secs(50)); // increase from default 5s

    let energy_model = ViennaRNA::default();
    let seq = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").unwrap();
    let db1 = DotBracketVec::try_from(".((((....)))).((((........))))...............").unwrap();
    let db2 = DotBracketVec::try_from(".((((....)))).((((.(....).))))...............").unwrap();
    let db3 = DotBracketVec::try_from(".((((....))))..(((........)))................").unwrap();
    let db4 = DotBracketVec::try_from(".((((....)))).((((.(.....)))))...............").unwrap();
    let db5 = DotBracketVec::try_from(".(((......))).((((........))))...............").unwrap();
    let db6 = DotBracketVec::try_from("..(((....)))..((((........))))...............").unwrap();
    let db7 = DotBracketVec::try_from(".(((......)))..(((........)))................").unwrap();
    let db8 = DotBracketVec::try_from(".(((.(...)))).((((........))))...............").unwrap();

    let macrostate = Macrostate::from_list(
        "bench_state",
        &seq,
        &[db1, db2, db3, db4, db5, db6, db7, db8], 
        &energy_model,
    );

    let rate_model = Metropolis::new(energy_model.temperature(), 1.0);
    group.bench_function("Neighborhood generation", |b| {
        b.iter(|| {
            let _ = ExitMacrostate::from((
                &macrostate, 
                &seq, 
                &energy_model, 
                &rate_model
            ));
        });
    });
}

criterion_group!(benches, exit_macrostates);
criterion_main!(benches);

