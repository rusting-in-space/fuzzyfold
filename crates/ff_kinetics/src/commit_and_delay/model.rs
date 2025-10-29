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
use crate::commit_and_delay::ExitMacrostateRegistry;

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

pub struct CommitAndDelay<'a, E: EnergyModel, R: RateModel> {
    exit_registry: Arc<ExitMacrostateRegistry<'a, E, R>>,
    trajectories: Array2<Option<ReactiveTrajectoryEnsemble>>,
}

impl<'a, E: EnergyModel, R: RateModel> From<Arc<ExitMacrostateRegistry<'a, E, R>>> 
for CommitAndDelay<'a, E, R> {
    fn from(exit_registry: Arc<ExitMacrostateRegistry<'a, E, R>>
    ) -> Self {
        let n = exit_registry.len();
        Self {
            exit_registry,
            trajectories: Array2::from_elem((n, n), None),
        }
    }
}

impl<'a, E: EnergyModel, R: RateModel> CommitAndDelay<'a, E, R> {

    pub fn simulate_from(&mut self, start_id: MacrostateID) {
        let sequence = self.exit_registry.parent_registry().sequence();
        let energy_model = self.exit_registry.parent_registry().energy_model();
        let rate_model = self.exit_registry.rate_model();

        let start_ms = self.exit_registry.exit_macrostates()
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
                let stop_id = self.exit_registry.parent_registry().classify(&stop_db);
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

    pub fn simulate_between(&mut self, start_id: MacrostateID, stop_id: MacrostateID) {
        let sequence = self.exit_registry.parent_registry().sequence();
        let energy_model = self.exit_registry.parent_registry().energy_model();
        let rate_model = self.exit_registry.rate_model();

        let start_ms = self.exit_registry.exit_macrostates()
            .get(start_id)
            .expect("invalid macrostate index");

        let start_db = start_ms.get_random_microstate().unwrap();
        let pairings = PairTable::try_from(&start_db).unwrap();
        let loops = LoopStructure::try_from((&sequence[..], &pairings, energy_model)).unwrap();
        let mut simulator = LoopStructureSSA::from((loops, rate_model));

        let mut mean_time = 0.0;
        let mut curr_db = start_db;
        let mut curr_id = start_id;
        let mut toggle = 0;
        simulator.simulate(
            &mut rng(), 
            f64::MAX,
            |t, _tinc, flux, ls| {
                let next_db = DotBracketVec::from(ls);
                let next_id = self.exit_registry.parent_registry().classify(&next_db);
                println!("current: {} {} {}", next_db, next_id, t);
                //TODO: think about this more..
                if next_id != toggle {
                    if next_id != 0 {
                        let traj = ReactiveMicroTrajectory {
                            i: curr_db.clone(),
                            j: next_db.clone(),
                            simu_time: t,
                            mean_time,
                        };
                        self.trajectories
                            .get_mut((curr_id, next_id))
                            .unwrap()
                            .get_or_insert_with(|| ReactiveTrajectoryEnsemble {
                                start: curr_id,
                                stop: next_id,
                                t_min: None,
                                t_max: None,
                                successes: Vec::new(),
                            })
                        .successes
                            .push(traj);
                        if next_id == stop_id {
                            return false;
                        } 
                        curr_db = next_db;
                        curr_id = next_id;
                        mean_time = 0.0;
                    }
                    toggle = next_id;

                }
                mean_time += 1.0/flux;
                true
            },
        );
    }

    pub fn gather_data(&mut self) {
        // let's make an executable. 
    }

    pub fn to_rate_matrix(&self) -> Array2<f64> {
        todo!("let's do only k_commit for now")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use ff_energy::ViennaRNA;
    use ff_energy::NucleotideVec;
    use crate::Metropolis;
    use crate::macrostates::MacrostateRegistry;

    fn test_ms1() -> std::io::Cursor<&'static [u8]> {
        Cursor::new(b">lmin=lm1_bh=4.0
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
        .(((.....((((.((((........))).).)))).....))).")
    }

    fn test_ms2() -> std::io::Cursor<&'static [u8]> {
        Cursor::new(b">lmin=lm2_bh=4.0
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....)))).((((.((((...))))..)))).........
        .((((....)))).((((.(((.....)))..)))).........
        ..(((....)))..((((.((((...))))..)))).........
        .((((....)))).(((..((((...))))...))).........
        .((((....))))..(((.((((...))))..)))..........
        .((((....)))).((((..(((...)))...)))).........
        .(((......))).((((.((((...))))..)))).........
        .((((....))))..(((.(((.....)))..)))..........
        ..(((....)))..((((.(((.....)))..)))).........
        .(((.(...)))).((((.((((...))))..)))).........
        .((((....))))..(((..(((...)))...)))..........
        .(((......))).((((.(((.....)))..)))).........
        ...((....))...((((.((((...))))..)))).........
        .(((......))).((((..(((...)))...)))).........
        .((((....)))).((((.((((...))).).)))).........
        .((((....)))).((((.((((...))))..))).)........
        ..(((....)))...(((.((((...))))..)))..........
        .(((......)))..(((.((((...))))..)))..........
        .((((....)))).((((...((...))....)))).........
        ..((......))..((((.((((...))))..)))).........
        .((((....)))).((((..((.....))...)))).........")
    }
 
    fn test_ms3() -> std::io::Cursor<&'static [u8]> {
        Cursor::new(b">lmin=lm3_bh=3.0
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....)))).((((........))))...............
        .((((....)))).((((.(....).))))...............
        .((((....))))..(((........)))................
        .((((....)))).((((.(.....)))))...............
        .(((......))).((((........))))...............
        ..(((....)))..((((........))))...............
        .(((......)))..(((........)))................
        .(((.(...)))).((((........))))...............")
    }

    #[test]
    fn test_commit_and_delay_minimal() {
        let energy_model = ViennaRNA::default();
        let seq = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").unwrap();
        let mut registry = MacrostateRegistry::from((&seq, &energy_model));

        registry.insert_from_reader(test_ms1(), "manual").unwrap();
        assert_eq!(registry.len(), 2);

       assert_eq!(registry.len(), 2);
        let rate_model = Metropolis::new(energy_model.temperature(), 1.0);
        let exitreg = ExitMacrostateRegistry::from((&registry, &rate_model));

        let mut cad = CommitAndDelay::from(Arc::new(exitreg));
        cad.simulate_from(1);
        assert_eq!(cad.trajectories.get((1, 1)).and_then(|opt| opt.as_ref()).unwrap().len(), 1);
        cad.simulate_from(1);
        assert_eq!(cad.trajectories.get((1, 1)).and_then(|opt| opt.as_ref()).unwrap().len(), 2);
        cad.simulate_from(1);
        assert_eq!(cad.trajectories.get((1, 1)).and_then(|opt| opt.as_ref()).unwrap().len(), 3);
    }

    #[test]
    fn test_commit_and_delay() {
        let energy_model = ViennaRNA::default();
        let seq = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").unwrap();
        let mut registry = MacrostateRegistry::from((&seq, &energy_model));

        registry.insert_from_reader(test_ms1(), "manual").unwrap();
        registry.insert_from_reader(test_ms2(), "manual").unwrap();
        registry.insert_from_reader(test_ms3(), "manual").unwrap();
        assert_eq!(registry.len(), 4);

        let rate_model = Metropolis::new(energy_model.temperature(), 1.0);
        let _exitreg = ExitMacrostateRegistry::from((&registry, &rate_model));

        //NOTE: too slow for a unittest at the moment.
        // let mut cad = CommitAndDelay::from(Arc::new(exitreg));
        // cad.simulate_between(3,1);
        // cad.simulate_between(1,3);
        // cad.simulate_between(2,1);
        // cad.simulate_between(1,2);
        // cad.simulate_between(3,2);
        // cad.simulate_between(2,3);
        // for row in cad.trajectories.rows() {
        //     let line = row
        //         .iter()
        //         .map(|el| el.as_ref().map_or("0".into(), |ens| ens.len().to_string()))
        //         .collect::<Vec<_>>()
        //         .join(" ");
        //     println!("{line}");
        // }
    }
}

