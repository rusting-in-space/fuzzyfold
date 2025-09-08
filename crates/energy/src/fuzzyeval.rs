use core::f64;
use std::i32;
use std::path::Path;


use crate::NearestNeighborLoop;
use crate::LoopDecomposition;
use crate::Base;
use crate::PairTypeRNA;
use crate::EnergyTables;
use crate::ParamError;
use crate::EnergyModel;
use crate::coaxial_stacking::CoaxialStackingDP;

fn rescale_energy_to_temp(enth: i32, en37: i32, temp_c: f64) -> i32 {
    let t_ref = 310.15; // 37 °C in Kelvin
    let t = temp_c + 273.15;
    let dtemp = t / t_ref;

    let enth = enth as f64;
    let en37 = en37 as f64;

    let dg = enth * (1.0 - dtemp) + en37 * dtemp;
    dg.round() as i32
}


#[derive(Debug)]
pub struct FuzzyEval {
    temperature: f64,
    min_hp_size: usize,

    _duplex_initiation_en37: i32,
    _duplex_initiation_enth: i32,
    terminal_ru_en37: i32,
    terminal_ru_enth: i32,
    _symmetry_en37: i32,
    _symmetry_enth: i32,

    coaxial_mm_discontious_en37: i32,
    coaxial_mm_discontious_enth: i32,
    coaxial_mm_wcf_bonus_en37: i32,
    coaxial_mm_wcf_bonus_enth: i32,
    coaxial_mm_gu_bonus_en37: i32,
    coaxial_mm_gu_bonus_enth: i32,

    energy_tables: EnergyTables,
}

impl FuzzyEval {
    pub fn from_parameter_file<P: AsRef<Path>>(path: P) -> Result<Self, ParamError> {
        let energy_tables = EnergyTables::from_parameter_file(path)?;
        Ok(FuzzyEval {
            temperature: 37.0,
            min_hp_size: 3,

            _duplex_initiation_en37: 410,
            _duplex_initiation_enth: 360,
            terminal_ru_en37: 50,
            terminal_ru_enth: 370,
            _symmetry_en37: 43,
            _symmetry_enth: 0,

            coaxial_mm_discontious_en37: -210,
            coaxial_mm_discontious_enth: -846, // +/- 2.75
            coaxial_mm_wcf_bonus_en37: -40,
            coaxial_mm_wcf_bonus_enth: -40,
            coaxial_mm_gu_bonus_en37: -20,
            coaxial_mm_gu_bonus_enth: -20,

            energy_tables,
        })
    }

    fn eval_hairpin_loop(&self, seq: &[Base]) -> i32 {
        let et = &self.energy_tables;
        let n = seq.len() - 2;
        if n < self.min_hp_size {
            panic!("Invalid hairpin size {n}");
        }
        if PairTypeRNA::from((seq[0], *seq.last().unwrap())) == PairTypeRNA::NN {
            panic!("Invalid closing pair {:?}", seq);
        }

        // Special hairpin energies
        if seq.len() <= 6 { 
            if let Some((en37, enth)) = et.hairpin_sequences.get(seq).copied() {
                if self.temperature == 37.0 { 
                    return en37;
                } else {
                    return rescale_energy_to_temp(enth, en37, self.temperature);
                }
            }
        }

        // Initiation terms
        let (mut en37, mut enth) = if n <= 30 {
            (et.hairpin[n].unwrap()
            ,et.hairpin_enthalpies[n].unwrap())
        } else {
            todo!("not implemented");
        };

        // NOTE: double check NN desrciption if this is indeed the correct way.
        if n == 3 && PairTypeRNA::from((seq[0], *seq.last().unwrap())).is_ru() {
            en37 += self.terminal_ru_en37;
            enth += self.terminal_ru_enth;
        } else if n > 3 {
            en37 += self.energy_tables.mismatch_hairpin
                [PairTypeRNA::from((seq[0], seq[n+1])) as usize]
                [seq[1] as usize][seq[n] as usize].unwrap();
            enth += self.energy_tables.mismatch_hairpin_enthalpies
                [PairTypeRNA::from((seq[0], seq[n+1])) as usize]
                [seq[1] as usize][seq[n] as usize].unwrap();
        }

        if self.temperature == 37.0 { 
            en37
        } else {
            rescale_energy_to_temp(enth, en37, self.temperature)
        }
    }

    /// g37 = g37 AU end-penalty + sum(37-stacking)
    fn eval_interior_loop(&self, fwdseq: &[Base], revseq: &[Base]) -> i32 {
        let outer = PairTypeRNA::from((*fwdseq.first().unwrap(), *revseq.last().unwrap()));
        let inner = PairTypeRNA::from((*revseq.first().unwrap(), *fwdseq.last().unwrap()));

        let (en37, enth) = match (fwdseq.len(), revseq.len()) {
            (2, 2) => (
                self.energy_tables.stack
                [outer as usize][inner as usize].expect("from file")
                ,
                self.energy_tables.stack_enthalpies
                [outer as usize][inner as usize].expect("from file")
            ),
            (3, 2) | (2, 3) => (//NOTE: SpecialC if C adjacent to paired C missing!
                self.energy_tables.bulge[1].unwrap() + 
                self.energy_tables.stack
                [outer as usize][inner as usize]
                .expect("from file")
                ,
                self.energy_tables.bulge_enthalpies[1].unwrap() +
                self.energy_tables.stack_enthalpies
                [outer as usize][inner as usize]
                .expect("from file")
            ),
            (3, 3) => (
                self.energy_tables.int11
                [outer as usize][inner as usize]
                [fwdseq[1] as usize][revseq[1] as usize].expect("from file")
                ,
                self.energy_tables.int11_enthalpies
                [outer as usize][inner as usize]
                [fwdseq[1] as usize][revseq[1] as usize].expect("from file")
            ),
            (4, 3) => (
                self.energy_tables.int21
                [outer as usize][inner as usize]
                [fwdseq[1] as usize][fwdseq[2] as usize]
                [revseq[1] as usize].expect("from file")
                ,
                self.energy_tables.int21_enthalpies
                [outer as usize][inner as usize]
                [fwdseq[1] as usize][fwdseq[2] as usize]
                [revseq[1] as usize].expect("from file")
            ),
            (3, 4) => (
                self.energy_tables.int21
                [inner as usize][outer as usize]
                [revseq[1] as usize][revseq[2] as usize]
                [fwdseq[1] as usize].expect("from file")
                ,
                self.energy_tables.int21_enthalpies
                [inner as usize][outer as usize]
                [revseq[1] as usize][revseq[2] as usize]
                [fwdseq[1] as usize].expect("from file")
            ),
            (4, 4) => (
                self.energy_tables.int22
                [outer as usize][inner as usize]
                [fwdseq[1] as usize][fwdseq[2] as usize]
                [revseq[1] as usize][revseq[2] as usize]
                .expect("from file")
                ,
                self.energy_tables.int22_enthalpies
                [outer as usize][inner as usize]
                [fwdseq[1] as usize][fwdseq[2] as usize]
                [revseq[1] as usize][revseq[2] as usize]
                .expect("from file")
            ),
            (l, 2) | (2, l) => { // General Bulge case
                let n = l - 2;
                let (pg1, th1) = if outer.is_ru() { (self.terminal_ru_en37, self.terminal_ru_enth) } else { (0, 0) };
                let (pg2, th2) = if inner.is_ru() { (self.terminal_ru_en37, self.terminal_ru_enth) } else { (0, 0) };
                if n <= 30 {(
                    self.energy_tables.bulge[n].unwrap() + pg1 + pg2,
                    self.energy_tables.bulge_enthalpies[n].unwrap() + th1 + th2
                )} else {
                    todo!("not implemented");
                }
            },
            (lfwd, lrev) => { 
                let n = (lfwd as isize - lrev as isize).abs() as usize;
                let (pg1, th1) = if outer.is_ru() { (self.terminal_ru_en37, self.terminal_ru_enth) } else { (0, 0) };
                let (pg2, th2) = if inner.is_ru() { (self.terminal_ru_en37, self.terminal_ru_enth) } else { (0, 0) };

                if n <= 30 {(
                    self.energy_tables.interior[n].unwrap() + pg1 + pg2,
                    self.energy_tables.interior_enthalpies[n].unwrap() + th1 + th2
                )} else {
                    todo!("not implemented");
                }
            }
        };

        if self.temperature == 37.0 { 
            en37
        } else {
            rescale_energy_to_temp(enth, en37, self.temperature)
        }
    }

    fn eval_multibranch_loop(&self, segments: &[&[Base]]) -> i32 {
        let n = segments.len(); 
        let mut en37 = 0;
        let mut enth = 0;
        let mut table = CoaxialStackingDP::new(n);
        for i in 0..n {
            let j = (i + 1) % n; 
            let pair = PairTypeRNA::from((*segments[i].last().unwrap(), segments[j][0]));
            if pair.is_ru() { 
                en37 += self.terminal_ru_en37;
                enth += self.terminal_ru_enth;
            }

            //println!("initE {i} | initS {i} {j} {:?}", pair);
            //let d5 = segments.get(i)
            //    .and_then(|seg| seg.len().checked_sub(2).and_then(|d| seg.get(d)));
            //let d3 = segments.get(j)
            //    .and_then(|seg| seg.len().checked_sub(2).and_then(|d| seg.get(1)));
 
            let d5 = d5_base(segments[i], false);
            let d3 = d3_base(segments[j], false);
            table.set_min_energy(i, d5, pair, d3, &self.energy_tables);

            let slen = segments[j].len() - 2;  
            if slen <= 1 {
                let j2 = (i + 2) % n; 
                let mm = if slen == 1 { segments[i].get(1) } else { None };
                let npair: PairTypeRNA = (*segments[j].last().unwrap(), segments[j2][0]).into();
                table.set_coax_energy(
                    i, j, d5, pair, mm, npair, d3, 
                    &self.energy_tables, 
                    self.coaxial_mm_discontious_en37);
 
            }
        }

        //let branches = segments.len() as i32;
        //let avg_asym = (2.0f64).min({
        //    let mut asy = Vec::new();
        //    for i in 0..segments.len() {
        //        let l = segments[i].len() as f64;
        //        let r = segments[i + 1 % segments.len()].len() as f64;
        //        asy.push((l - r).abs());
        //    }
        //    assert!(asy.len() > 1);
        //    asy.iter().sum::<f64>() / asy.len() as f64
        //});

        if self.temperature == 37.0 { 
            en37 + table.compute_multibranch_mfe(segments)
        } else {
            rescale_energy_to_temp(enth, en37, self.temperature)
        }
    }

    fn eval_exterior_loop(&self, segments: &[&[Base]]) -> i32 {
        let n = segments.len() - 1; 
        let mut en37 = 0;
        let mut enth = 0;
        let mut table = CoaxialStackingDP::new(n);
        for i in 0..n {
            let j = i + 1; 
            let pair = PairTypeRNA::from((*segments[i].last().unwrap(), segments[j][0]));
            if pair.is_ru() { 
                en37 += self.terminal_ru_en37;
                enth += self.terminal_ru_enth;
            }

            let d5 = d5_base(segments[i], i == 0);
            let d3 = d3_base(segments[j], j == n);
            table.set_min_energy(i, d5, pair, d3, &self.energy_tables);

            if j >= n {
                continue;
            }

            let slen = segments[j].len() - 2;  
            if slen <= 1 {
                let j2 = i + 2; 
                let mm = if slen == 1 { segments[i].get(1) } else { None };
                let npair = (*segments[j].last().unwrap(), segments[j2][0]).into();
                table.set_coax_energy(
                    i, j, d5, pair, mm, npair, d3, 
                    &self.energy_tables, 
                    self.coaxial_mm_discontious_en37);
            }
        }

        if self.temperature == 37.0 { 
            en37 + table.compute_exterior_mfe(segments)
        } else {
            //TODO
            rescale_energy_to_temp(enth, en37, self.temperature)
        }
    }

}

impl EnergyModel for FuzzyEval {

    fn can_pair(&self, b1: Base, b2: Base) -> bool {
        PairTypeRNA::from((b1, b2)).can_pair()
    }

    fn min_hairpin_size(&self) -> usize { self.min_hp_size }

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

    fn energy_of_loop(&self, sequence: &[Base], nn_loop: &NearestNeighborLoop) -> i32 {

        match nn_loop {
            NearestNeighborLoop::Hairpin { closing: (i, j) } => {
                self.eval_hairpin_loop(&sequence[*i..=*j])
            }
            NearestNeighborLoop::Interior { closing: (i, j), inner: (k, l) } => {
                let left = &sequence[*i..=*k];
                let right = &sequence[*l..=*j];
                self.eval_interior_loop(left, right)
            }
            NearestNeighborLoop::Multibranch { closing: (i, j), branches } => {
                let mut slices: Vec<&[Base]> = Vec::new();
                let mut start = *i;
                for &(k, l) in branches {
                    slices.push(&sequence[start..=k]);
                    start = l;
                }
                slices.push(&sequence[start..=*j]);
                self.eval_multibranch_loop(&slices)
            }
            NearestNeighborLoop::Exterior { branches } => {
                let mut slices: Vec<&[Base]> = Vec::new();
                let mut start = 0;
                for &(k, l) in branches {
                    slices.push(&sequence[start..=k]);
                    start = l;
                }
                slices.push(&sequence[start..]);
                self.eval_exterior_loop(&slices)
            }
        }
    }
}

pub fn d5_base(seg: &[Base], exterior: bool) -> Option<&Base> {
    if exterior && seg.len() >= 2 {
        seg.get(seg.len() - 2)
    } else if !exterior && seg.len() >= 3 {
        seg.get(seg.len() - 2)
    } else {
        None
    }
}

pub fn d3_base(seg: &[Base], exterior: bool) -> Option<&Base> {
    if exterior && seg.len() >= 2 {
        seg.get(1)
    } else if !exterior && seg.len() >= 3 {
        seg.get(1)
    } else {
        None
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::basify;

    #[test]
    fn test_ff_hairpin_evaluation() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/params/rna_turner2004.par");
        let model = FuzzyEval::from_parameter_file(path).unwrap();

        assert_eq!(model.eval_hairpin_loop(&basify("GAAAC")), 540);
    }

    #[test]
    fn test_ff_exterior_evaluation() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/params/rna_turner2004.par");
        let model = FuzzyEval::from_parameter_file(path).unwrap();

        let seg1 = basify("UUG");
        let seg2 = basify("CUG");
        let seg3 = basify("CG");
        let seg4 = basify("CUG");

        let binding: Vec<&[Base]> = vec![&seg1, &seg2, &seg3, &seg4];

        assert_eq!(model.eval_exterior_loop(&binding), -450); //TODO

    }

    #[test]
    fn test_ff_multibranch_evaluation() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/params/rna_turner2004.par");
        let model = FuzzyEval::from_parameter_file(path).unwrap();

        let seg1 = basify("UUG");
        let seg2 = basify("CUG");
        let seg3 = basify("CG");
        let seg4 = basify("CUG");

        let binding: Vec<&[Base]> = vec![&seg1, &seg2, &seg3, &seg4];

        assert_eq!(model.eval_multibranch_loop(&binding), -610); //TODO

    }


}

 
