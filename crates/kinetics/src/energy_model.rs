// Intuition:
//let model = ViennaRNA::from_parameter_file("turner2004.par");
//let energy = model.energy_of_structure(&sequence, &pairings);
//let energy = model.energy_of_pair(&sequence, &pairings, i, j);

use crate::utils::Base;
use crate::NearestNeighborLoop;
use crate::LoopDecomposition;
use structure::PairTable;

#[derive(Debug)]
pub struct ViennaRNA {
    dangles: u8,
    temperature: i32, // Kelvin?
}

impl ViennaRNA {
    pub fn new() -> Self {
        ViennaRNA {
            dangles: 0,
            temperature: 37,
        }
    }

    fn hairpin(&self, _seq: &[Base]) -> i32 {
        //assert!(can_pair(*seq[0], seq.last().unwrap(), true));
        -3
    }

    fn interior(&self, _fwseq: &[Base], _revseq: &[Base]) -> i32 {
       -2 
    }

    fn multibranch(&self, segments: &[&[Base]]) -> i32 {
        -1 * (segments.len() as i32)
    }

    fn exterior(&self, _segments: &[&[Base]]) -> i32 {
        0 
    }
}

pub trait EnergyModel {
    fn can_pair(&self, b1: Base, b2: Base) -> bool;

    #[inline] fn min_hp_size(&self) -> usize { 3 }

    fn energy_of_structure<T: LoopDecomposition>(&self, 
        sequence: &[Base], 
        structure: &T
    ) -> i32;

    /// The free energy difference between forming the pair vs breaking it.
    //TODO: may not be particularly efficient, and depends on PairTable!
    fn energy_of_pair(&self, 
        sequence: &[Base], 
        structure: &PairTable,
        i: usize,
        j: usize,
    ) -> i32;

    fn energy_of_loop(&self, 
        sequence: &[Base], 
        nn_loop: &NearestNeighborLoop
    ) -> i32;
}

impl EnergyModel for ViennaRNA {

    fn can_pair(&self, b1: Base, b2: Base) -> bool {
        matches!((b1, b2),
        (Base::A, Base::U) | (Base::U, Base::A) |
        (Base::G, Base::C) | (Base::C, Base::G) |
        (Base::G, Base::U) | (Base::U, Base::G))
    }

    fn energy_of_structure<T: LoopDecomposition>(&self, 
        sequence: &[Base], 
        structure: &T
    ) -> i32  {
        let mut total = 0;
        structure.for_each_loop(|l| {
            total += self.energy_of_loop(sequence, l);
        });
        total
    }

    fn energy_of_pair(&self, 
        sequence: &[Base], 
        structure: &PairTable, 
        i: usize,
        j: usize
    ) -> i32 {
        let inner = structure.loop_enclosed_by(Some((i, j)));
        let epair = structure.get_enclosing_pair(i, j);

        if let (Some(pi), Some(pj)) = (structure[i], structure[j]) {
            assert!(j == pi && i == pj);
            let outer = structure.loop_enclosed_by(epair);
            let mut pt = structure.clone();
            pt[i] = None;
            pt[j] = None;
            let combo = pt.loop_enclosed_by(epair);

            let en_paired = self.energy_of_loop(sequence, &inner) + self.energy_of_loop(sequence, &outer);
            let en_absent = self.energy_of_loop(sequence, &combo);
            return en_absent - en_paired
        } else {
            assert!(structure[i] == structure[j]);
            let combo = structure.loop_enclosed_by(epair);
            let mut pt = structure.clone();
            pt[i] = Some(j);
            pt[j] = Some(i);
            let outer = pt.loop_enclosed_by(epair);
            let en_paired = self.energy_of_loop(sequence, &inner) + self.energy_of_loop(sequence, &outer);
            let en_absent = self.energy_of_loop(sequence, &combo);
            return en_paired - en_absent 
        }
    }

    fn energy_of_loop(&self, sequence: &[Base], nn_loop: &NearestNeighborLoop) -> i32 {

        match nn_loop {
            NearestNeighborLoop::Hairpin { closing: (i, j) } => {
                let loop_seq = &sequence[i + 1..*j];
                self.hairpin(loop_seq)
            }
            NearestNeighborLoop::Interior { closing: (i, j), inner: (k, l) } => {
                let left = &sequence[*i..=*k];
                let right = &sequence[*l..=*j];
                self.interior(left, right)
            }
            NearestNeighborLoop::Multibranch { closing: (i, j), branches } => {
                let mut slices: Vec<&[Base]> = Vec::new();
                let mut start = *i;
                for &(k, l) in branches {
                    slices.push(&sequence[start..=k]);
                    start = l;
                }
                slices.push(&sequence[start..=*j]);
                self.multibranch(&slices)
            }
            NearestNeighborLoop::Exterior { branches } => {
                let mut slices: Vec<&[Base]> = Vec::new();
                let mut start = 0;
                for &(k, l) in branches {
                    slices.push(&sequence[start..=k]);
                    start = l;
                }
                slices.push(&sequence[start..]);
                self.exterior(&slices)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::basify;

    #[test]
    fn test_basic_stuff() {
        let sequence = "GCAUACGAUCA";
        let struct_1 = ".(....)....";

        let model = ViennaRNA::new();
        let sequence = basify(sequence);
        let struct_1 = PairTable::try_from(struct_1).expect("valid");

        let energy_1 = model.energy_of_structure(&sequence, &struct_1);
        println!("Total energy: {}", energy_1);

        let energy = model.energy_of_pair(&sequence, &struct_1, 1, 6);
        println!("Pair energy: {}", energy);
    }


    #[test]
    fn test_fake_model() {
        let sequence = "GCAUACGAUCA";
        let pairings = ".((.(.)()))";

        let model = ViennaRNA::new();
        let sequence = basify(sequence);
        let pairings = PairTable::try_from(pairings).expect("valid");

        let energy = model.energy_of_structure(&sequence, &pairings);
        println!("Total energy: {}", energy);

        let energy = model.energy_of_pair(&sequence, &pairings, 2, 9);
        println!("Pair energy: {}", energy);
    }
}

 
