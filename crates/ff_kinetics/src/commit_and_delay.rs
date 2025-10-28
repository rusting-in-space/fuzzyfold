use std::sync::Arc;
use std::convert::TryFrom;
use rand::rng;
use ndarray::Array2;
use ff_structure::PairTable;
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;

use crate::RateModel;
use crate::LoopStructure;
use crate::LoopStructureSSA;
use crate::macrostates::ExitMacrostateRegistry;

type MacrostateID = usize;

#[derive(Debug, Clone)]
struct ReactiveMicroTrajectory {
    i: DotBracketVec,
    j: DotBracketVec,
    simu_time: f64,
    mean_time: f64,
}

#[derive(Debug, Clone)]
struct ReactiveTrajectoryEnsemble {
    start: MacrostateID,
    stop: MacrostateID,
    t_min: Option<f64>,
    t_max: Option<f64>,
    successes: Vec<ReactiveMicroTrajectory>,
}

impl ReactiveTrajectoryEnsemble {
    pub fn sort_trajectories(&mut self) {
        self.successes.sort_by(|a, b| a.simu_time.partial_cmp(&b.simu_time).unwrap());
    }
    
    pub fn split_ensemble(&self, num_splits: usize) -> Vec<Self> {
        if self.successes.is_empty() {
            return vec![];
        }

        //TODO: test with num_splits = 0
        let chunk_size = self.successes.len().div_ceil(num_splits);
        self.successes
            .chunks(chunk_size)
            .map(|chunk| Self {
                start: self.start,
                stop: self.stop,
                t_min: chunk.first().map(|t| t.simu_time),
                t_max: chunk.last().map(|t| t.simu_time),
                successes: chunk.to_vec(),
            })
            .collect()
    }

    pub fn len(&self) -> usize {
        self.successes.len()
    }
}

struct CommitAndDelay<'a, E: EnergyModel, R: RateModel> {
    pub registry: Arc<ExitMacrostateRegistry<'a, E, R>>,
    pub trajectories: Array2<Option<ReactiveTrajectoryEnsemble>>,
}

impl<'a, E: EnergyModel, R: RateModel> From<Arc<ExitMacrostateRegistry<'a, E, R>>> 
for CommitAndDelay<'a, E, R> {
    fn from(registry: Arc<ExitMacrostateRegistry<'a, E, R>>
    ) -> Self {
        let n = registry.registry.len();
        Self {
            registry,
            trajectories: Array2::from_elem((n, n), None),
        }
    }
}


impl<'a, E: EnergyModel, R: RateModel> CommitAndDelay<'a, E, R> {

    pub fn simulate_from(&mut self, start_id: MacrostateID) {
        let sequence = self.registry.registry.sequence;
        let energy_model = self.registry.registry.energy_model;
        let rate_model = self.registry.rate_model;

        let start_ms = self.registry.exit_macrostates
            .get(start_id)
            .expect("invalid macrostate index");

        let start_db = start_ms.get_random_microstate().unwrap();
        let pairings = PairTable::try_from(&start_db).unwrap();
        let loops = LoopStructure::try_from((&sequence[..], &pairings, energy_model)).unwrap();
        let mut simulator = LoopStructureSSA::from((loops, rate_model));

        let mut mean_time = 0.0;
        simulator.simulate(
            &mut rng(), 
            f64::MAX,
            |t, _tinc, flux, ls| {
                let stop_db = DotBracketVec::from(ls);
                let stop_id = self.registry.registry.classify(&stop_db);
                //println!("current: {} {} {}", stop_db, stop_id, t);
                if stop_id != 0usize {
                    let traj = ReactiveMicroTrajectory {
                        i: start_db.clone(),
                        j: stop_db,
                        simu_time: t,
                        mean_time,
                    };
                    self.trajectories
                        .get_mut((start_id, stop_id))
                        .unwrap()
                        .get_or_insert_with(|| ReactiveTrajectoryEnsemble {
                            start: start_id,
                            stop: stop_id,
                            t_min: None,
                            t_max: None,
                            successes: Vec::new(),
                        })
                        .successes
                        .push(traj);
                    return false;
                }
                mean_time += 1.0/flux;
                true
            },
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ff_energy::ViennaRNA;
    use ff_energy::NucleotideVec;
    use crate::macrostates::Macrostate;
    use crate::Metropolis;
    use crate::macrostates::MacrostateRegistry;

    #[test]
    fn test_commit_and_delay() {
        /*
        >lmin=lm1_bh=4.0
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....((((.((((........))))..))))....)))).
        ..(((....((((.((((........))))..))))....)))..
        .((((....((((.((((.(....).))))..))))....)))).
        .(.((....((((.((((........))))..))))....)).).
        .((((....((((.((((.(.....)))))..))))....)))).
        .((.(....((((.((((........))))..))))....).)).
        .((((....(((..((((........))))...)))....)))).
        .((((....((((..(((........)))...))))....)))).
        .((((....((((.(((..........)))..))))....)))).
        .((((.....(((.((((........))))..))).....)))).
        .(((.....((((.((((........))))..)))).....))).
        .(((.....((((..(((........)))...)))).....))).
        ...((....((((.((((........))))..))))....))...
        .((((....((((.((((........))).).))))....)))).
        ..((.....((((.((((........))))..)))).....))..
        .((((....((((.((((........)))).).)))....)))).
        .(((.....((((.((((.(.....)))))..)))).....))).
        .((......((((.((((........))))..))))......)).
        .(((......(((.((((........))))..)))......))).
        .(((.....((((.((((........))).).)))).....))).
        */

        let seq = NucleotideVec::from_lossy("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC");
        let energy_model = ViennaRNA::default();
        let rate_model = Metropolis::new(energy_model.temperature(), 1.0);

        let mut registry = MacrostateRegistry::from((&seq, &energy_model));

        // Build a test macrostate with a few structures
        let s1 = DotBracketVec::try_from(".((((....((((.((((........))))..))))....)))).").unwrap(); 
        let s2 = DotBracketVec::try_from("..(((....((((.((((........))))..))))....)))..").unwrap();
        let s3 = DotBracketVec::try_from(".((((....((((.((((.(....).))))..))))....)))).").unwrap();
        let s4 = DotBracketVec::try_from(".(.((....((((.((((........))))..))))....)).).").unwrap();
        let s5 = DotBracketVec::try_from(".((((....((((.((((.(.....)))))..))))....)))).").unwrap();
        let s6 = DotBracketVec::try_from(".((.(....((((.((((........))))..))))....).)).").unwrap();
        //let s7 = DotBracketVec::try_from(".((((....(((..((((........))))...)))....)))).").unwrap();
        //let s8 = DotBracketVec::try_from(".((((....((((..(((........)))...))))....)))).").unwrap();
        //let s9 = DotBracketVec::try_from(".((((....((((.(((..........)))..))))....)))).").unwrap();
        //let s10 = DotBracketVec::try_from(".((((.....(((.((((........))))..))).....)))).").unwrap();
        //let s11 = DotBracketVec::try_from(".(((.....((((.((((........))))..)))).....))).").unwrap();
        //let s12 = DotBracketVec::try_from(".(((.....((((..(((........)))...)))).....))).").unwrap();
        //let s13 = DotBracketVec::try_from("...((....((((.((((........))))..))))....))...").unwrap();
        //let s14 = DotBracketVec::try_from(".((((....((((.((((........))).).))))....)))).").unwrap();
        //let s15 = DotBracketVec::try_from("..((.....((((.((((........))))..)))).....))..").unwrap();
        //let s16 = DotBracketVec::try_from(".((((....((((.((((........)))).).)))....)))).").unwrap();
        //let s17 = DotBracketVec::try_from(".(((.....((((.((((.(.....)))))..)))).....))).").unwrap();
        //let s18 = DotBracketVec::try_from(".((......((((.((((........))))..))))......)).").unwrap();
        //let s19 = DotBracketVec::try_from(".(((......(((.((((........))))..)))......))).").unwrap();
        //let s20 = DotBracketVec::try_from(".(((.....((((.((((........))).).)))).....))).").unwrap();
 

        let macrostate = Macrostate::from_list(
            "test",
            &seq,
            &[s1.clone(), s2.clone(), s3.clone(),
              s4.clone(), s5.clone(), s6.clone()],
            &energy_model,
        );
        registry.macrostates.push(macrostate);

        let exitreg = ExitMacrostateRegistry::from((&registry, &rate_model));

        let mut cad = CommitAndDelay::from(Arc::new(exitreg));
        cad.simulate_from(1);
        assert_eq!(cad.trajectories.get((1, 1)).and_then(|opt| opt.as_ref()).unwrap().len(), 1);
        cad.simulate_from(1);
        assert_eq!(cad.trajectories.get((1, 1)).and_then(|opt| opt.as_ref()).unwrap().len(), 2);
        cad.simulate_from(1);
        assert_eq!(cad.trajectories.get((1, 1)).and_then(|opt| opt.as_ref()).unwrap().len(), 3);
    }
}

