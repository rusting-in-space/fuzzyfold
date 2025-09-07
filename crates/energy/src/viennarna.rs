use core::f64;
use std::i32;
use std::path::Path;
use structure::PairTable;

use crate::NearestNeighborLoop;
use crate::LoopDecomposition;
use crate::Base;
use crate::PairTypeRNA;
use crate::EnergyTables;
use crate::ParamError;
use crate::EnergyModel;

fn rescale_energy_to_temp(enth: i32, en37: i32, temp_c: f64) -> i32 {
    let t_ref = 310.15; // 37 °C in Kelvin
    let t = temp_c + 273.15;
    let dtemp = t / t_ref;

    let enth = enth as f64;
    let en37 = en37 as f64;

    let dg = enth * (1.0 - dtemp) + en37 * dtemp;
    dg as i32 // No rounding for ViennaRNA compatibility
}


/// The default ViennaRNA-v2.6 energy model.
///
/// The current implementation allows different
/// parameter files, but only the standard
/// settings (dangle=2, tetraloops, etc.)
///
/// Only single-stranded folding is supported.
///
#[derive(Debug)]
pub struct ViennaRNA {
    temperature: f64,
    min_hp_size: usize,

    lxc37: f64, /* ViennaRNA parameter for logarithmic loop energy extrapolation */

    // ML params section -- hardcoded.
    ml_base_en37: i32,
    ml_base_enth: i32,
    ml_closing_en37: i32,
    ml_closing_enth: i32,
    ml_intern_en37: i32,
    ml_intern_enth: i32,

    ninio_en37: i32,
    ninio_enth: i32,
    max_ninio: i32,

    _duplex_initiation_en37: i32,
    _duplex_initiation_enth: i32,
    terminal_ru_en37: i32,
    terminal_ru_enth: i32,

    energy_tables: EnergyTables,
}

impl ViennaRNA {

    pub fn default() -> Self {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/params/rna_turner2004.par");
        ViennaRNA::from_parameter_file(path)
            .expect("Built-in Turner 2004 parameter file must be valid")
    }

    pub fn from_parameter_file<P: AsRef<Path>>(path: P) -> Result<Self, ParamError> {
        let energy_tables = EnergyTables::from_parameter_file(path)?;
        Ok(ViennaRNA {
            temperature: 37.0,
            min_hp_size: 3,

            lxc37: 107.856, //TODO
                            
            // ML params section -- hardcoded.
            ml_base_en37: 0,
            ml_base_enth: 0,
            ml_closing_en37: 930,
            ml_closing_enth: 3000,
            ml_intern_en37: -90,
            ml_intern_enth: -220,

            // NINIO params section -- hardcoded.
            ninio_en37: 60,
            ninio_enth: 320,
            max_ninio: 300,

            // Misc params section -- hardcoded.
            _duplex_initiation_en37: 410,
            _duplex_initiation_enth: 360,
            terminal_ru_en37: 50,
            terminal_ru_enth: 370,

            energy_tables,
        })
    }

    fn hairpin(&self, seq: &[Base]) -> Result<i32, ParamError> {
        let n = seq.len() - 2;
        if n < self.min_hp_size {
            return Err(ParamError::InvalidHairpinSize(n));
        }

        // Special hairpin energies
        if seq.len() <= 6 { 
            if let Some((en37, enth)) = self.energy_tables.hairpin_sequences.get(seq).copied() {
                if self.temperature == 37.0 { 
                    return Ok(en37);
                } else {
                    return Ok(rescale_energy_to_temp(enth, en37, self.temperature))
                }
            }
        }

        // Initiation terms
        let (mut en37, mut enth) = if n <= 30 {
            (self.energy_tables.hairpin[n].ok_or(ParamError::MissingValue("hairpin", n))?,
            self.energy_tables.hairpin_enthalpies[n].ok_or(ParamError::MissingValue("hairpin", n))?)
        } else {
            (self.energy_tables.hairpin[30].ok_or(ParamError::MissingValue("hairpin", 30))?
                + (self.lxc37 * (n as f64).ln() / 30.) as i32,
             self.energy_tables.hairpin_enthalpies[30].ok_or(ParamError::MissingValue("hairpin", 30))?
                + (self.lxc37 * (n as f64).ln() / 30.) as i32)
        };

        // NOTE: double check NN desrciption if this is indeed the correct way.
        if n == 3 && matches!((seq[0], *seq.last().unwrap()).into()
            , PairTypeRNA::GU | PairTypeRNA::UG 
            | PairTypeRNA::AU | PairTypeRNA::UA
            | PairTypeRNA::NN
        ) { 
            en37 += self.terminal_ru_en37;
            enth += self.terminal_ru_enth;
        } else if n > 3 {
            en37 += self.energy_tables.mismatch_hairpin[
                PairTypeRNA::from((seq[0], seq[n+1])) as usize][ 
                    seq[1] as usize][
                    seq[n] as usize
                    ].ok_or(ParamError::MissingValue("mismatch_hairpin", n))?;

            enth += self.energy_tables.mismatch_hairpin_enthalpies[
                PairTypeRNA::from((seq[0], seq[n+1])) as usize][ 
                    seq[1] as usize][
                    seq[n] as usize
                    ].ok_or(ParamError::MissingValue("mismatch_hairpin", n))?;
        }

        if self.temperature == 37.0 { 
            Ok(en37)
        } else {
            Ok(rescale_energy_to_temp(enth, en37, self.temperature))
        }
    }

    /// g37 = g37 AU end-penalty + sum(37-stacking)
    fn interior(&self, fwdseq: &[Base], revseq: &[Base]) -> i32 {
        let outer = PairTypeRNA::from((*fwdseq.first().unwrap(), *revseq.last().unwrap()));
        let inner = PairTypeRNA::from((*revseq.first().unwrap(), *fwdseq.last().unwrap()));

        let is_ru_end = |pt| matches!(pt
            , PairTypeRNA::GU | PairTypeRNA::UG 
            | PairTypeRNA::AU | PairTypeRNA::UA | PairTypeRNA::NN);
        

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
                let (pg1, th1) = if is_ru_end(outer) { (self.terminal_ru_en37, self.terminal_ru_enth) } else { (0, 0) };
                let (pg2, th2) = if is_ru_end(inner) { (self.terminal_ru_en37, self.terminal_ru_enth) } else { (0, 0) };
                if n <= 30 {(
                    self.energy_tables.bulge[n].unwrap() + pg1 + pg2,
                    self.energy_tables.bulge_enthalpies[n].unwrap() + th1 + th2
                )} else {(
                    self.energy_tables.bulge[30].unwrap() + pg1 + pg2 +
                    (self.lxc37 * (n as f64).ln() / 30.) as i32,
                    self.energy_tables.bulge_enthalpies[30].unwrap() + th1 + th2 +
                    (self.lxc37 * (n as f64).ln() / 30.) as i32
                )}
            },
            (lfwd, lrev) => { 
                let n = (lfwd as isize - lrev as isize).abs() as usize;
                let (pg1, th1) = if is_ru_end(outer) { (self.terminal_ru_en37, self.terminal_ru_enth) } else { (0, 0) };
                let (pg2, th2) = if is_ru_end(inner) { (self.terminal_ru_en37, self.terminal_ru_enth) } else { (0, 0) };

                if n <= 30 {(
                    self.energy_tables.interior[n].unwrap() + pg1 + pg2,
                    self.energy_tables.interior_enthalpies[n].unwrap() + th1 + th2
                )} else {(
                    self.energy_tables.interior[30].unwrap() + pg1 + pg2 +
                    (self.lxc37 * (n as f64).ln() / 30.) as i32,
                    self.energy_tables.interior_enthalpies[30].unwrap() + th1 + th2 +
                    (self.lxc37 * (n as f64).ln() / 30.) as i32
                )}
            }
        };

        if self.temperature == 37.0 { 
            en37
        } else {
            rescale_energy_to_temp(enth, en37, self.temperature)
        }
    }

    fn multibranch(&self, segments: &[&[Base]]) -> i32 {
        let is_ru_end = |pt| matches!(pt
            , PairTypeRNA::GU | PairTypeRNA::UG 
            | PairTypeRNA::AU | PairTypeRNA::UA | PairTypeRNA::NN);

        let n = segments.len(); 

        let mut en37 = 0;
        let mut enth = 0;
 
        for i in 0..n {
            let j = (i + 1) % n; 
            let pair = PairTypeRNA::from((*segments[i].last().unwrap(), segments[j][0]));
            if is_ru_end(pair) { 
                en37 += self.terminal_ru_en37;
                enth += self.terminal_ru_enth;
            }
            let d5 = segments.get(i)
                .and_then(|seg| seg.len().checked_sub(2).and_then(|d| seg.get(d)));
            let d3 = segments.get(j)
                .and_then(|seg| seg.len().checked_sub(2).and_then(|d| seg.get(1)));

            //NOTE: This does not take the minimum over all options, it always
            // prefers terminal mismatch over single dangling.
            let (g, h) = match (d5, d3) { 
                (Some(&b5), Some(&b3)) => 
                    (self.energy_tables.mismatch_exterior[pair as usize][b5 as usize][b3 as usize].unwrap(),
                     self.energy_tables.mismatch_exterior_enthalpies[pair as usize][b5 as usize][b3 as usize].unwrap()),
                (Some(&b5), None) => 
                    (self.energy_tables.dangle5[pair as usize][b5 as usize].unwrap(),
                     self.energy_tables.dangle5_enthalpies[pair as usize][b5 as usize].unwrap()),
                (None, Some(&b3)) => 
                    (self.energy_tables.dangle3[pair as usize][b3 as usize].unwrap(),
                     self.energy_tables.dangle3_enthalpies[pair as usize][b3 as usize].unwrap()),
                _ => (0, 0),
            };
            en37 += g;
            enth += h;
        }
 
        en37 += self.ml_base_en37 
             + self.ml_closing_en37 
             + self.ml_intern_en37 * n as i32;
        enth += self.ml_base_enth 
             + self.ml_closing_enth 
             + self.ml_intern_enth * n as i32;

        if self.temperature == 37.0 { 
            en37
        } else {
            rescale_energy_to_temp(enth, en37, self.temperature)
        }
    }

    fn exterior(&self, segments: &[&[Base]]) -> i32 {
        let is_ru_end = |pt| matches!(pt
            , PairTypeRNA::GU | PairTypeRNA::UG 
            | PairTypeRNA::AU | PairTypeRNA::UA | PairTypeRNA::NN);
 
        let mut en37 = 0;
        let mut enth = 0;
        let n = segments.len() - 1; 
        for i in 0..n {
            let j = i + 1; 
            
            let pair = (*segments[i].last().unwrap(), segments[j][0]).into();
            if is_ru_end(pair) { 
                en37 += self.terminal_ru_en37;
                enth += self.terminal_ru_enth;
            }

            let d5 = segments.get(i)
                .and_then(|seg| seg.len().checked_sub(2).and_then(|d| seg.get(d)));
            let d3 = segments.get(j)
                .and_then(|seg| seg.len().checked_sub(2).and_then(|d| seg.get(1)));

            //NOTE: This does not take the minimum over all options, it always
            // prefers terminal mismatch over single dangling.
            let (g, h) = match (d5, d3) { 
                (Some(&b5), Some(&b3)) => 
                    (self.energy_tables.mismatch_exterior[pair as usize][b5 as usize][b3 as usize].unwrap(),
                     self.energy_tables.mismatch_exterior_enthalpies[pair as usize][b5 as usize][b3 as usize].unwrap()),
                (Some(&b5), None) => 
                    (self.energy_tables.dangle5[pair as usize][b5 as usize].unwrap(),
                     self.energy_tables.dangle5_enthalpies[pair as usize][b5 as usize].unwrap()),
                (None, Some(&b3)) => 
                    (self.energy_tables.dangle3[pair as usize][b3 as usize].unwrap(),
                     self.energy_tables.dangle3_enthalpies[pair as usize][b3 as usize].unwrap()),
                _ => (0, 0),
            };
            en37 += g;
            enth += h;
        }
        if self.temperature == 37.0 { 
            en37
        } else {
            rescale_energy_to_temp(enth, en37, self.temperature)
        }
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

}

impl EnergyModel for ViennaRNA {

    fn can_pair(&self, b1: Base, b2: Base) -> bool {
        matches!((b1, b2),
        (Base::A, Base::U) | (Base::U, Base::A) |
        (Base::G, Base::C) | (Base::C, Base::G) |
        (Base::G, Base::U) | (Base::U, Base::G))
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
    use crate::basify;

    #[test]
    fn test_vrna_hairpin_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.hairpin(&basify("GAAAC")).unwrap(), 540);
        assert_eq!(model.hairpin(&basify("CCGAGG")).unwrap(), 350);
        assert_eq!(model.hairpin(&basify("CCAAGG")).unwrap(), 330);
        assert_eq!(model.hairpin(&basify("CAAGG")).unwrap(), 540);
        assert_eq!(model.hairpin(&basify("CAAAG")).unwrap(), 540);
        assert_eq!(model.hairpin(&basify("AAAAU")).unwrap(), 590);
        assert_eq!(model.hairpin(&basify("GAAAU")).unwrap(), 590);
        assert_eq!(model.hairpin(&basify("CAAAAG")).unwrap(), 410);
        assert_eq!(model.hairpin(&basify("ACCCU")).unwrap(), 590);
        assert_eq!(model.hairpin(&basify("GCCCCC")).unwrap(), 490);

        assert_eq!(model.hairpin(&basify("GAAAG")).unwrap(), 590);
        assert_eq!(model.hairpin(&basify("CAAAC")).unwrap(), 590);

        assert_eq!(model.hairpin(&basify("AAAAAU")).unwrap(), 530);
        assert_eq!(model.hairpin(&basify("GAAAAU")).unwrap(), 580);
        assert_eq!(model.hairpin(&basify("ACCCCU")).unwrap(), 540);
        assert_eq!(model.hairpin(&basify("ACCCCCU")).unwrap(), 550);
        assert_eq!(model.hairpin(&basify("AAAAAAU")).unwrap(), 540);
        assert_eq!(model.hairpin(&basify("AAAAAAAU")).unwrap(), 510);
        assert_eq!(model.hairpin(&basify("AAAAAAAAAAU")).unwrap(), 610);
        assert_eq!(model.hairpin(&basify("AAAAAAAAAAA")).unwrap(), 660);
    }

    #[test]
    fn test_vrna_stacking_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.interior(&basify("CG"), &basify("CG")), -240);
        assert_eq!(model.interior(&basify("AC"), &basify("GU")), -220);
        assert_eq!(model.interior(&basify("GU"), &basify("AC")), -220);
    }

    #[test]
    fn test_vrna_int11_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.interior(&basify("CCG"), &basify("CGG")), 50);
        assert_eq!(model.interior(&basify("CAG"), &basify("CAG")), 90);
        assert_eq!(model.interior(&basify("ACU"), &basify("AAU")), 190);
        assert_eq!(model.interior(&basify("GCU"), &basify("AUC")), 120);
        assert_eq!(model.interior(&basify("GCU"), &basify("AGC")), 120);
    }

    #[test]
    fn test_vrna_int21_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.interior(&basify("CACG"), &basify("CGG")), 110);
        assert_eq!(model.interior(&basify("CAAG"), &basify("CAG")), 230);
        assert_eq!(model.interior(&basify("AACU"), &basify("AAU")), 370);
        assert_eq!(model.interior(&basify("GACU"), &basify("AUC")), 300);
        assert_eq!(model.interior(&basify("GACU"), &basify("AGC")), 300);

        assert_eq!(model.interior(&basify("CGG"), &basify("CACG")), 110);
        assert_eq!(model.interior(&basify("CAG"), &basify("CAAG")), 230);
        assert_eq!(model.interior(&basify("AAU"), &basify("AACU")), 370);
        assert_eq!(model.interior(&basify("AUC"), &basify("GACU")), 300);
        assert_eq!(model.interior(&basify("AGC"), &basify("GACU")), 300);
    }

    #[test]
    fn test_vrna_bulge_1_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.interior(&basify("CAG"), &basify("CG")), 140);
        assert_eq!(model.interior(&basify("AAU"), &basify("AU")), 270);
        assert_eq!(model.interior(&basify("GAU"), &basify("AC")), 160);

        assert_eq!(model.interior(&basify("CCG"), &basify("CG")), 140);
        assert_eq!(model.interior(&basify("ACU"), &basify("AU")), 270);
        assert_eq!(model.interior(&basify("GCU"), &basify("AC")), 160);

        assert_eq!(model.interior(&basify("CG"), &basify("CAG")), 140);
        assert_eq!(model.interior(&basify("AU"), &basify("AAU")), 270);
        assert_eq!(model.interior(&basify("AC"), &basify("GAU")), 160);
                                                                
        assert_eq!(model.interior(&basify("CG"), &basify("CCG")), 140);
        assert_eq!(model.interior(&basify("AU"), &basify("ACU")), 270);
        assert_eq!(model.interior(&basify("AC"), &basify("GCU")), 160);
    }

    #[test]
    fn test_vrna_bulge_2_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.interior(&basify("CAAG"), &basify("CG")), 280);
        assert_eq!(model.interior(&basify("AAAU"), &basify("AU")), 380);
        assert_eq!(model.interior(&basify("GAAU"), &basify("AC")), 330);

        assert_eq!(model.interior(&basify("CCAG"), &basify("CG")), 280);
        assert_eq!(model.interior(&basify("ACAU"), &basify("AU")), 380);
        assert_eq!(model.interior(&basify("GCAU"), &basify("AC")), 330);

        assert_eq!(model.interior(&basify("CG"), &basify("CAAG")), 280);
        assert_eq!(model.interior(&basify("AU"), &basify("AAAU")), 380);
        assert_eq!(model.interior(&basify("AC"), &basify("GAAU")), 330);

        assert_eq!(model.interior(&basify("CG"), &basify("CCAG")), 280);
        assert_eq!(model.interior(&basify("AU"), &basify("ACAU")), 380);
        assert_eq!(model.interior(&basify("AC"), &basify("GCAU")), 330);
    }

    #[test]
    fn test_vrna_bulge_n_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.interior(&basify("CAAAAAAG"), &basify("CG")), 440);
        assert_eq!(model.interior(&basify("CAAAAAAAAG"), &basify("CG")), 470);
    }

    #[test]
    fn test_vrna_multibranch() {
        let model = ViennaRNA::default();

        let seg1 = &basify("AAAA");
        let seg2 = &basify("AAA");
        let seg3 = &basify("AAAAA");
        let energy = model.multibranch(&[seg1, seg2, seg3]);
        assert_eq!(energy, 720);

        let seg1 = &basify("GAAC");
        let seg2 = &basify("GAC");
        let seg3 = &basify("GAAAC");
        let energy = model.multibranch(&[seg1, seg2, seg3]);
        assert_eq!(energy, 330);

        let seg1 = &basify("GAAC");
        let seg2 = &basify("GAC");
        let seg3 = &basify("GAAAAAAAAAAAAAAAAAAC");
        let energy = model.multibranch(&[seg1, seg2, seg3]);
        assert_eq!(energy, 330);

    }

    #[test]
    fn test_vrna_exterior_single_branch() {
        let model = ViennaRNA::default();

        let seg1 = &basify("AUG");
        let seg2 = &basify("CUG");
        let energy = model.exterior(&[seg1, seg2]);
        assert_eq!(energy, -120);

        let seg1 = &basify("UG");
        let seg2 = &basify("CU");
        let energy = model.exterior(&[seg1, seg2]);
        assert_eq!(energy, -120); 

        let seg1 = &basify("G");
        let seg2 = &basify("CU");
        let energy = model.exterior(&[seg1, seg2]);
        assert_eq!(energy, -120);
 
        let seg1 = &basify("UG");
        let seg2 = &basify("C");
        let energy = model.exterior(&[seg1, seg2]);
        assert_eq!(energy, 0); 
    }

    #[test]
    fn test_vrna_exterior_single_branch_t25() {
        let mut model = ViennaRNA::default();
        model.temperature = 25.0;

        let seg1 = &basify("AUG");
        let seg2 = &basify("CUG");
        let energy = model.exterior(&[seg1, seg2]);
        assert_eq!(energy, -148);

        let seg1 = &basify("UG");
        let seg2 = &basify("CU");
        let energy = model.exterior(&[seg1, seg2]);
        assert_eq!(energy, -148);

        let seg1 = &basify("G");
        let seg2 = &basify("CU");
        let energy = model.exterior(&[seg1, seg2]);
        assert_eq!(energy, -144);
 
        let seg1 = &basify("UG");
        let seg2 = &basify("C");
        let energy = model.exterior(&[seg1, seg2]);
        assert_eq!(energy, -1); 
 
        let seg1 = &basify("A");
        let seg2 = &basify("U");
        let energy = model.exterior(&[seg1, seg2]);
        assert_eq!(energy, 62); 
 
    }

    #[test]
    fn test_vrna_exterior_two_branches() {
        let model = ViennaRNA::default();

        let seg1 = &basify("AUG");
        let seg2 = &basify("CUG");
        let seg3 = &basify("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3]);
        assert_eq!(energy, -240);

        let seg1 = &basify("AUG");
        let seg2 = &basify("CUUG");
        let seg3 = &basify("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3]);
        assert_eq!(energy, -240);

        let seg1 = &basify("AUG");
        let seg2 = &basify("CUUG");
        let seg3 = &basify("C");
        let energy = model.exterior(&[seg1, seg2, seg3]);
        assert_eq!(energy, -120);

        let seg1 = &basify("AUG");
        let seg2 = &basify("CG");
        let seg3 = &basify("CU");
        let energy = model.exterior(&[seg1, seg2, seg3]);
        assert_eq!(energy, -290);

        let seg1 = &basify("ACA");
        let seg2 = &basify("UGG");
        let seg3 = &basify("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3]);
        assert_eq!(energy, -130);
    }

    #[test]
    fn test_vrna_exterior_three_branches() {
        let model = ViennaRNA::default();

        let seg1 = &basify("AUG");
        let seg2 = &basify("CUG");
        let seg3 = &basify("CUG");
        let seg4 = &basify("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3, seg4]);
        assert_eq!(energy, -360);

        let seg1 = &basify("AUG");
        let seg2 = &basify("CUG");
        let seg3 = &basify("UUG");
        let seg4 = &basify("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3, seg4]);
        assert_eq!(energy, -240);
    }

 
    #[test]
    fn test_evaluations() {
        let model = ViennaRNA::default();

        let seq = "GAAAAC";
        let dbr = "(....)";
        let e37 = 450;
        assert_eq!(model.energy_of_structure(&basify(seq), &PairTable::try_from(dbr).expect("valid")), e37);

        let seq = "GAAAAG";
        let dbr = "(....)";
        let e37 = 630;
        assert_eq!(model.energy_of_structure(&basify(seq), &PairTable::try_from(dbr).expect("valid")), e37);

        let seq = "ACGUUAAAGACGU";
        let dbr = "(((((...)))))";
        let e37 = -170;
        assert_eq!(model.energy_of_structure(&basify(seq), &PairTable::try_from(dbr).expect("valid")), e37);
    }
}

 
