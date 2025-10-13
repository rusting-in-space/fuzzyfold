
use crate::NAIDX;
use crate::PairTable;

pub trait Neighborhood {
    fn add_neighbors(&self) -> impl Iterator<Item = (NAIDX, NAIDX)> + '_;
    fn del_neighbors(&self) -> impl Iterator<Item = (NAIDX, NAIDX)> + '_;
}

impl Neighborhood for PairTable {
    /// All *deletion* moves: every existing base pair (i,j)
   }
    fn del_neighbors(&self) -> impl Iterator<Item = (NAIDX, NAIDX)> + '_ {
        self.0.iter().enumerate()
            .filter_map(|(i, &partner)| match partner {
                Some(j) if j > i as NAIDX => Some((i as NAIDX, j)),
                _ => None,
            })
    }


    /// All *insertion* moves: possible new base pairs (i,j)
    fn add_neighbors(&self) -> impl Iterator<Item = (NAIDX, NAIDX)> + '_ {
        (0..self.len())
            .filter(|&i| self[i].is_none())
            .flat_map(move |i| {
                (i + 1..self.len())
                    .scan(None, move |_, j| {
                        if let Some(pj) = self[j] {
                            let pj = pj as usize;
                            if pj < i {
                                // outer closing pair -> next loop
                                return None;
                            }
                            if pj > j {
                                // skip paired segment (neighboring helix)
                                return Some(j + (pj - j) + 1);
                            }
                        }
                        Some(j)
                    })
                    .filter(move |&j| self[j].is_none())
                    .map(move |j| (i as NAIDX, j as NAIDX))
            })
    }
}


