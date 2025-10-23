use ff_structure::NAIDX;
use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;

use crate::RateModel;
use crate::LoopStructure;

pub trait ApplyMove {
    fn apply_move(&mut self, mv: Move);
    fn undo_move(&mut self, mv: Move) {
        self.apply_move(mv.inverse());
    }
}

impl ApplyMove for DotBracketVec {
    fn apply_move(&mut self, mv: Move) {
        match mv {
            Move::Add { i, j } => {
                self[i as usize] = DotBracket::Open;
                self[j as usize] = DotBracket::Close;
            }
            Move::Del { i, j } => {
                self[i as usize] = DotBracket::Unpaired;
                self[j as usize] = DotBracket::Unpaired;
            }
        }
    }
}

impl<'a, M: EnergyModel> ApplyMove for LoopStructure<'a, M> { 
    fn apply_move(&mut self, mv: Move) {
        match mv {
            Move::Add { i, j } => {
                self.apply_add_move(i, j);
            }
            Move::Del { i, j } => {
                self.apply_del_move(i, j);
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Move {
    Add {
        i: NAIDX,
        j: NAIDX,
    },
    Del {
        i: NAIDX,
        j: NAIDX,
    },
}

impl Move {
    pub fn inverse(self) -> Self {
        match self {
            Move::Add { i, j } => Move::Del { i, j },
            Move::Del { i, j } => Move::Add { i, j },
        }
    }
}


#[derive(Debug, Clone, PartialEq)]
pub enum Reaction {
    Add {
        i: NAIDX,
        j: NAIDX,
        delta_e: i32,
        log_rate: f64,
    },
    Del {
        i: NAIDX,
        j: NAIDX,
        delta_e: i32,
        log_rate: f64,
    },
}

impl Reaction {
    pub fn new_add<K: RateModel>(model: &K, 
        i: NAIDX, j: NAIDX, delta_e: i32
) -> Self {
        let rate = model.log_rate(delta_e);
        Reaction::Add { i, j, delta_e, log_rate: rate }
    }

    pub fn new_del<K: RateModel>(model: &K, 
        i: NAIDX, j: NAIDX, delta_e: i32) -> Self {
        let rate = model.log_rate(delta_e);
        Reaction::Del { i, j, delta_e, log_rate: rate }
    }

    pub fn ij(&self) -> (NAIDX, NAIDX) {
        match self {
            Reaction::Add { i, j, .. } => (*i, *j),
            Reaction::Del { i, j, .. } => (*i, *j),
        }
    }

    pub fn log_rate(&self) -> f64 {
        match self {
            Reaction::Add { log_rate, .. } => *log_rate,
            Reaction::Del { log_rate, .. } => *log_rate,
        }
    }

    pub fn delta_e(&self) -> i32 {
        match self {
            Reaction::Add { delta_e, .. } => *delta_e,
            Reaction::Del { delta_e, .. } => *delta_e,
        }
    }

}

