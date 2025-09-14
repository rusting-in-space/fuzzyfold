
use std::f64;

use rand::Rng;

use energy::NucleotideVec;
use energy::ViennaRNA;
use energy::EnergyModel;
use structure::PairTable;
use kinetics::LoopStructure;

/// One possible move in the system
#[derive(Debug, Clone)]
pub enum Reaction {
    Add {
        i: usize,
        j: usize,
        delta_e: i32, // energy difference in model units
        rate: f64,
    },
    Del {
        i: usize,
        j: usize,
        delta_e: i32,
        rate: f64,
    },
}

impl Reaction {
    pub fn rate(&self) -> f64 {
        match self {
            Reaction::Add { rate, .. } => *rate,
            Reaction::Del { rate, .. } => *rate,
        }
    }
}

/// Convert ΔE into a rate using Metropolis rule, k0 = 1
fn metropolis_rate(delta_e: i32, rt: f64, k0: f64) -> f64 {
    let delta = delta_e as f64 / 100.0; // assume your model returns energies in 0.01 kcal/mol
    if delta <= 0.0 {
        k0
    } else {
        k0 * (-delta / rt).exp()
    }
}

/// Collect all reactions (adds + dels) from a LoopStructure
fn collect_reactions<M: EnergyModel>(
    ls: &LoopStructure<M>,
    rt: f64,
    k0: f64,
) -> Vec<Reaction> {
    let mut reactions = Vec::new();

    for (i, j, delta) in ls.get_add_neighbors() {
        let rate = metropolis_rate(delta, rt, k0);
        reactions.push(Reaction::Add { i, j, delta_e: delta, rate });
    }

    for (i, j, delta) in ls.get_del_neighbors() {
        let rate = metropolis_rate(delta, rt, k0);
        reactions.push(Reaction::Del { i, j, delta_e: delta, rate });
    }

    reactions
}

/// Run Gillespie SSA until t_max, target structure reached, or no more moves
pub fn run_ssa<M: EnergyModel>(
    mut ls: LoopStructure<M>,
    rt: f64,
    k0: f64,
    t_max: f64,
    _target_structures: Option<Vec<PairTable>>,
) -> (f64, LoopStructure<M>) {
    let mut rng = rand::rng();
    let mut t = 0.0;

    println!("{:>15.9} {}", t, ls);

    loop {
        // 1. Collect reactions
        let reactions = collect_reactions(&ls, rt, k0);
        if reactions.is_empty() {
            panic!("An emptpy set of reactions is no good.");
        }

        // 2. Total rate
        let a0: f64 = reactions.iter().map(|r| r.rate()).sum();
        if a0 <= 0.0 {
            panic!("An flux a <= 0 is no good.");
        }

        // 3. Sample tau
        let tau = -(rng.random::<f64>().ln()) / a0;
        if t + tau > t_max {
            t = t_max;
            println!("{:>15.9} {}", t, ls);
            break;
        }
        t += tau;

        // 4. Sample which reaction
        let mut r = rng.random::<f64>() * a0;
        let mut chosen: Option<Reaction> = None;
        for rxn in reactions {
            r -= rxn.rate();
            if r <= 0.0 {
                chosen = Some(rxn);
                break;
            }
        }
        let chosen = chosen.expect("no reaction chosen");

        // 5. Apply chosen reaction
        match chosen {
            Reaction::Add { i, j, .. } => { let _ = ls.apply_add_move(i, j); },
            Reaction::Del { i, j, .. } => { let _ = ls.apply_del_move(i, j); },
        }

        println!("{:>15.9} {}", t, ls);
    }

    (t, ls)
}

/// Example main driving the SSA
fn main() {
    // Example setup: sequence, structure, and model
    let sequence = NucleotideVec::from_lossy("GCCCCGGUCA");
    let structure =      PairTable::try_from("...........").unwrap();
    let model = ViennaRNA::default();
    let rt = 0.001987 * (273.15 + 37.0); // ~0.616 kcal/mol
    let k0 = 1.;

    let ls = LoopStructure::try_from((&sequence[..], &structure, &model)).unwrap();
    let _ = run_ssa(ls, rt, k0, 1e7, None);

}

