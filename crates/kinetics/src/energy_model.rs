
use std::path::Path;

use crate::NearestNeighborLoop;
use crate::LoopDecomposition;
use structure::PairTable;

use crate::energy_tables::Base;
use crate::energy_tables::pair_type;
use crate::energy_tables::EnergyTables;
use crate::energy_tables::ParamError;

#[derive(Debug)]
pub struct ViennaRNA {
    temperature: f64,
    min_hp_size: usize,
    lxc37: f64, /* ViennaRNA parameter for logarithmic loop energy extrapolation */
    energy_tables: EnergyTables,

    dangles: u8,
}

impl ViennaRNA {
    pub fn from_parameter_file<P: AsRef<Path>>(path: P) -> Result<Self, ParamError> {
        let energy_tables = EnergyTables::from_parameter_file(path)?;
        Ok(ViennaRNA {
            temperature: 37.0,
            min_hp_size: 3,
            lxc37: 107.856, //TODO
            energy_tables,

            dangles: 0,
        })
    }

    fn hairpin(&self, seq: &[Base]) -> Result<i32, ParamError> {
        let n = seq.len() - 2;
        if n < self.min_hp_size {
            return Err(ParamError::InvalidHairpinSize(n));
        }

        if seq.len() <= 6 {
            if let Some((g37, _h)) = self.energy_tables.hairpin_sequences.get(seq).copied() {
                return Ok(g37);
            }
        }

        let mut energy = if n <= 30 {
            self.energy_tables.hairpin[n].ok_or(ParamError::MissingValue("hairpin", n))?
        } else {
            self.energy_tables.hairpin[30].ok_or(ParamError::MissingValue("hairpin", 30))?
                + (self.lxc37 * (n as f64).ln() / 30.) as i32
        };

        energy += self.energy_tables.mismatch_hairpin[
            pair_type(seq[0], seq[n+1]) as usize][ 
                seq[1] as usize][
                seq[n] as usize
                ].ok_or(ParamError::MissingValue("mismatch_hairpin", n))?;

        Ok(energy)
    }

    fn interior(&self, fwdseq: &[Base], revseq: &[Base]) -> i32 {
        if fwdseq.len() == 2 && revseq.len() == 2 {
           return self.energy_tables.stack[
               pair_type(fwdseq[0], revseq[1]) as usize][
               pair_type(revseq[0], fwdseq[0]) as usize].expect("from file"); 
        } else if fwdseq.len() == 3 && revseq.len() == 2 {
            return 0; // TODO: Bulge
        } else if fwdseq.len() == 2 && revseq.len() == 3 {
            return 0; // TODO: Bulge
        } else if fwdseq.len() == 3 && revseq.len() == 3 {
            // 1-1 interior loop
            return self.energy_tables.int11
                [pair_type(fwdseq[0], revseq[2]) as usize]
                [pair_type(revseq[0], fwdseq[2]) as usize]
                [fwdseq[1] as usize]
                [revseq[1] as usize]
                .unwrap(); 
        } else if fwdseq.len() == 4 && revseq.len() == 3 {
           return self.energy_tables.int21
               [pair_type(fwdseq[0], revseq[2]) as usize]
               [pair_type(revseq[0], fwdseq[3]) as usize]
               [fwdseq[1] as usize]
               [fwdseq[2] as usize]
               [revseq[1] as usize]
               .unwrap(); 
        } else if fwdseq.len() == 3 && revseq.len() == 4 {
           return self.energy_tables.int21
               [pair_type(revseq[0], fwdseq[2]) as usize]
               [pair_type(fwdseq[0], revseq[3]) as usize]
               [revseq[1] as usize]
               [revseq[2] as usize]
               [fwdseq[1] as usize]
               .unwrap(); 
        } else if fwdseq.len() == 4 && revseq.len() == 4 {
           return self.energy_tables.int22
               [pair_type(revseq[0], fwdseq[1]) as usize]
               [pair_type(fwdseq[0], revseq[2]) as usize]
               [fwdseq[1] as usize]
               [fwdseq[2] as usize]
               [revseq[1] as usize]
               [revseq[2] as usize]
               .unwrap(); 
        } 


        let n1 = fwdseq.len() - 2;
        let n2 = revseq.len() - 2;
        let a = -200.0; // 1.0 kcal/mol
        let b =   30.0; // 0.3 kcal/mol
        let c =   20.0; // 0.2 kcal/mol
        let energy = a + b * (n1 + n2) as f64 + c * (n1 as f64 - n2 as f64).abs();
        energy.round() as i32
    }

    fn multibranch(&self, segments: &[&[Base]]) -> i32 {
        let branches = segments.len();
        let unpaired: usize = segments.iter().map(|s| s.len() - 2).sum();
        let a = 340.0; // 3.4 kcal/mol = 34 dcal
        let b =  40.0; // 0.4 kcal/mol = 4 dcal
        let c =  10.0; // 0.1 kcal/mol = 1 dcal
        let energy = a + b * (branches as f64) + c * (unpaired as f64);
        energy.round() as i32
    }

    fn exterior(&self, segments: &[&[Base]]) -> i32 {
        if self.dangles == 0 || segments.len() == 0 {
            0 
        } else {
            segments.iter().map(|s| if s.len() - 2 > 0 { -50i32 } else { 0i32 }).sum::<i32>()
        } 
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
                self.hairpin(&sequence[*i..=*j]).unwrap()
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
    use crate::energy_tables::basify;

    #[test]
    fn test_basic_stuff() {
        let sequence = "GCAUACGAUCA";
        let struct_1 = ".(....)....";

        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/rna_turner2004.par");
        let model = ViennaRNA::from_parameter_file(path).unwrap();

        let sequence = basify(sequence);
        let struct_1 = PairTable::try_from(struct_1).expect("valid");

        let energy_1 = model.energy_of_structure(&sequence, &struct_1);
        println!("Total energy: {}", energy_1);

        let energy = model.energy_of_pair(&sequence, &struct_1, 1, 6);
        println!("Pair energy: {}", energy);
    }


    #[test]
    fn test_fake_model() {
        let sequence = "GCAUCCCCGAAAAUUG";
        let pairings = ".((.(...)(...)))";

        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/rna_turner2004.par");
        let model = ViennaRNA::from_parameter_file(path).unwrap();

        let sequence = basify(sequence);
        let pairings = PairTable::try_from(pairings).expect("valid");

        let energy = model.energy_of_structure(&sequence, &pairings);
        println!("Total energy: {}", energy);

        let energy = model.energy_of_pair(&sequence, &pairings, 2, 14);
        println!("Pair energy: {}", energy);
    }
}

 
