
use rand::rng;

use energy::NucleotideVec;
use energy::ViennaRNA;
use structure::PairTable;
use kinetics::LoopStructure;
use kinetics::LoopStructureSSA;
use kinetics::Metropolis;

/// Example main driving the SSA
fn main() {
    // Example setup: sequence, structure, and model
    let sequence = NucleotideVec::from_lossy("GCCCCGGUCA");
    let structure =      PairTable::try_from("...........").unwrap();
    let model = ViennaRNA::default();
    let rt = 0.001987 * (273.15 + 37.0); // ~0.616 kcal/mol
    let k0 = 1.;

    // Build LoopStructure from sequence + initial structure
    let loops = LoopStructure::try_from((&sequence[..], &structure, &model)).unwrap();
    let ratemodel = Metropolis::new(rt, k0);
    let mut simulator = LoopStructureSSA::from((loops, &ratemodel));
    simulator.simulate(
        &mut rng(), 
        0., 
        1e7, 
        |t, ls| {
            println!("{:15.9} {}", t, ls);
        }
    );

}

