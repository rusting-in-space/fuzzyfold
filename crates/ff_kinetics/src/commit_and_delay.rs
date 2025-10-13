use crate::macrostates::MacrostateRegistry;

type MicrostateID = usize;
type MacrostateID = usize;

struct ReactiveTrajectory {
    i: MicrostateID,
    j: MicrostateID,
    simu_time: f64,
    mean_time: f64,
}

struct ReactiveEnsemble {
    start: MacrostateID,
    stop: MacrostateID,
    successes: Vec<ReactiveTrajectory>,
}

impl ReactiveEnsemble {
    pub fn sort_tajectories(&mut self) {
    }

    pub fn get_subchannels(&self, num_splits: usize) -> Vec<ReactiveSubEnsemble> {
        let subchannels: Vec<ReactiveSubEnsemble> = Vec::with_capacity(num_splits);
        subchannels
    }
}

struct ReactiveSubEnsemble {
    start: MacrostateID,
    stop: MacrostateID,
    t_min: f64,
    t_max: f64,
    successes: Vec<ReactiveTrajectory>,
}

struct CommitAndDelay {
    macrostates: MacrostateRegistry,
}

impl CommitAndDelay {
    pub fn initialize(macrostates: MacrostateRegistry) {

    }

    fn trajectories() {

    }

    pub fn simulate_from(macrostate: MacrostateID, num: usize) {


    }

    pub fn simulate_all(num: usize) {

    }

}
