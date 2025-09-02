
use std::fmt;

use crate::energy_tables::EnergyTables;
use crate::ViennaRNA;
use crate::EnergyModel;
use crate::energy_tables::Base;
use crate::energy_tables::PairTypeRNA;
use crate::energy_tables::pair_type;
use crate::energy_tables::rev_pair_type;

struct DPTable {
    n: usize,
    data: Vec<i32>,
}

impl DPTable {
    fn new(n: usize) -> Self {
        DPTable {
            n,
            data: vec![i32::MAX / 2; n * n * 4], // d=4 for (a,b) states
        }
    }

    fn idx(&self, i: usize, j: usize, d5: usize, d3: usize) -> usize {
        (((i * self.n) + j) << 2) + ((d5 << 1) | d3)
    }

    fn get(&self, i: usize, j: usize, d5: usize, d3: usize) -> i32 {
        self.data[self.idx(i, j, d5, d3)]
    }

    fn set(&mut self, i: usize, j: usize, d5: usize, d3: usize, val: i32) {
        let idx = self.idx(i, j, d5, d3);
        self.data[idx] = val;
    }
}

fn set_min_energy(
    min_energy: &mut DPTable, 
    i: usize, 
    d5: Option<&Base>, 
    pi: PairTypeRNA,
    d3: Option<&Base>,
    etables: &EnergyTables,
) {
    min_energy.set(i, i, 0, 0, 0);
    min_energy.set(i, i, 1, 0, 
        d5.map_or(i32::MAX/2, |&b| {
            etables.dangle5[pi as usize][b as usize].unwrap()
        })
    );
    min_energy.set(i, i, 0, 1, 
        d3.map_or(i32::MAX/2, |&b| {
            etables.dangle3[pi as usize][b as usize].unwrap()
        })
    );
    min_energy.set(i, i, 1, 1, 
        match (d5, d3) { 
            (Some(&b5), Some(&b3)) => etables
                .mismatch_exterior[pi as usize]
                [b5 as usize][b3 as usize].unwrap(),
            _ => i32::MAX/2,
        }
    );
}
 
fn set_coax_energy(
    coax_energy: &mut DPTable, 
    i: usize, 
    d5: Option<&Base>, 
    pi: PairTypeRNA,
    mm: Option<&Base>, 
    pj: PairTypeRNA,
    d3: Option<&Base>,
    model: &ViennaRNA,
) {
    if let Some(&m) = mm {
        coax_energy.set(i, i + 1, 1, 0, 
            d5.map_or(i32::MAX/2, |&b| {
                let bonus = if matches!(pair_type(b, m), 
                    PairTypeRNA::GU | PairTypeRNA::UG) {
                    -20
                } else if model.can_pair(b, m) {
                    -40
                } else {
                    0
                };
                model.energy_tables.mismatch_exterior
                    [pi as usize][b as usize][m as usize].unwrap()
                    + model.g37_coaxial_mm_discontious
                        + bonus
            }));
        coax_energy.set(i, i + 1, 0, 1, 
            d3.map_or(i32::MAX/2, |&b| {
                let bonus = if matches!(pair_type(m, b), 
                    PairTypeRNA::GU | PairTypeRNA::UG) {
                    -20
                } else if model.can_pair(m, b) {
                    -40
                } else {
                    0
                };
                model.energy_tables.mismatch_exterior
                    [pi as usize][m as usize][b as usize].unwrap()
                    + model.g37_coaxial_mm_discontious
                        + bonus
            }));
    } else {
        coax_energy.set(i, i + 1, 0, 0, {
            let pj = rev_pair_type(&pj);
            model.energy_tables.stack[pi as usize][pj as usize].unwrap()
        });
    }
}

fn get_mfe(
    min_energy: &mut DPTable, 
    coax_energy: &mut DPTable, 
    ilen: usize, 
    i: usize, 
    jlen: usize, 
    j: usize
) {
    if jlen >= 2 {
        for &(d5, d3) in &[(0usize, 0usize), (0, 1), (1, 0), (1, 1)] {
            let best = (min_energy.get(i, j - 1, d5, 1) + min_energy.get(j, j, 1, d3))
                   .min(min_energy.get(i, j - 1, d5, 1) + min_energy.get(j, j, 0, d3))
                   .min(min_energy.get(i, j - 1, d5, 0) + min_energy.get(j, j, 1, d3))
                   .min(min_energy.get(i, j - 1, d5, 0) + min_energy.get(j, j, 0, d3));
            min_energy.set(i, j, d5, d3, best);
        }
    } else if jlen == 1 {
        for &(d5, d3) in &[(0usize, 0usize), (0, 1), (1, 0), (1, 1)] {
            let mut best = (min_energy.get(i, j - 1, d5, 1) + min_energy.get(j, j, 0, d3))
                   .min(min_energy.get(i, j - 1, d5, 0) + min_energy.get(j, j, 1, d3))
                   .min(min_energy.get(i, j - 1, d5, 0) + min_energy.get(j, j, 0, d3));
            if j == i + 1 {
                best = best.min(coax_energy.get(i, j, d5, d3));
            } else if ilen >= 2 {
                best = best
                    .min(min_energy.get(i, j - 2, d5, 1) + coax_energy.get(j-1, j, 1, d3))
                    .min(min_energy.get(i, j - 2, d5, 1) + coax_energy.get(j-1, j, 0, d3))
                    .min(min_energy.get(i, j - 2, d5, 0) + coax_energy.get(j-1, j, 1, d3))
                    .min(min_energy.get(i, j - 2, d5, 0) + coax_energy.get(j-1, j, 0, d3));
            } else if ilen == 1 {
                best = best
                    .min(min_energy.get(i, j - 2, d5, 1) + coax_energy.get(j-1, j, 0, d3))
                    .min(min_energy.get(i, j - 2, d5, 0) + coax_energy.get(j-1, j, 1, d3))
                    .min(min_energy.get(i, j - 2, d5, 0) + coax_energy.get(j-1, j, 0, d3));
            } else if ilen == 0 {
                best = best.min(min_energy.get(i, j - 2, d5, 0) + coax_energy.get(j-1, j, 0, d3));
            }
            min_energy.set(i, j, d5, d3, best);
        }
    } else if jlen == 0 {
        for &(d5, d3) in &[(0usize, 0usize), (0, 1), (1, 0), (1, 1)] {
            let mut best = min_energy.get(i, j - 1, d5, 0) + min_energy.get(j, j, 0, d3);
            if j == i + 1 {
                best = best.min(coax_energy.get(i, j, d5, d3));
            } else {
                if ilen >= 2 {
                    best = best
                        .min(min_energy.get(i, j - 2, d5, 1) + coax_energy.get(j-1, j, 1, d3))
                        .min(min_energy.get(i, j - 2, d5, 1) + coax_energy.get(j-1, j, 0, d3))
                        .min(min_energy.get(i, j - 2, d5, 0) + coax_energy.get(j-1, j, 1, d3))
                        .min(min_energy.get(i, j - 2, d5, 0) + coax_energy.get(j-1, j, 0, d3));
                } else if ilen == 1 {
                    best = best
                        .min(min_energy.get(i, j - 2, d5, 1) + coax_energy.get(j-1, j, 0, d3))
                        .min(min_energy.get(i, j - 2, d5, 0) + coax_energy.get(j-1, j, 1, d3))
                        .min(min_energy.get(i, j - 2, d5, 0) + coax_energy.get(j-1, j, 0, d3));
                } else if ilen == 0 {
                    best = best.min(min_energy.get(i, j - 2, d5, 0) + coax_energy.get(j-1, j, 0, d3));
                }
            }
            min_energy.set(i, j, d5, d3, best);
        }
    }
}

impl fmt::Debug for DPTable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "DPTable (n = {})", self.n)?;
        for i in 0..self.n {
            for j in i..self.n {
                // show all 4 context states
                let vals: Vec<String> = (0..=1)
                    .flat_map(|d5| (0..=1).map(move |d3| (d5, d3)))
                    .map(|(d5, d3)| {
                        let v = self.get(i, j, d5, d3);
                        if v >= i32::MAX / 4 {
                            format!("({},{}): ∞", d5, d3)
                        } else {
                            format!("({},{}): {}", d5, d3, v)
                        }
                    })
                    .collect();
                writeln!(f, "[{},{}] -> {}", i, j, vals.join(", "))?;
            }
        }
        Ok(())
    }
}

/// Main DP state
pub struct CoaxialStackingDP<'a> {
    segments: &'a [&'a [Base]],
    branching_pairs: Vec<PairTypeRNA>,
    min_energy: DPTable,  // E^{ab}[i..j]
    coax_energy: DPTable, // S^{ab}[i..j]
    model: &'a ViennaRNA,
}

impl<'a> CoaxialStackingDP<'a> {
 
    pub fn from_multibranch_loop(segments: &'a [&'a [Base]], model: &'a ViennaRNA) -> Self {
        let n = segments.len();
        let etables = &model.energy_tables;
        let mut branching_pairs = Vec::with_capacity(n);
        let mut min_energy = DPTable::new(n);
        let mut coax_energy = DPTable::new(n);
        for i in 0..n {
            let j = (i + 1) % n; 
            let pair = pair_type(*segments[i].last().unwrap(), segments[j][0]);
            println!("me {i} {j} {:?}", pair);
            branching_pairs.push(pair);
            let d5 = d5_base(segments[i], false);
            let d3 = d3_base(segments[j], false);
            set_min_energy(&mut min_energy, i, d5, pair, d3, etables);

            let slen = segments[j].len() - 2;  
            if slen <= 1 {
                let j2 = (i + 2) % n; 
                let mm = if slen == 1 { segments[i].get(1) } else { None };
                let npair = pair_type(*segments[j].last().unwrap(), segments[j2][0]);
                println!("cs {j} {j2} {:?}", npair);
                set_coax_energy(&mut coax_energy, i, d5, pair, mm, npair, d3, model);
            }
            // We now have E_i and S_ij.
            // So we could calculate E_0..i!
        }

        for i in 0..=1 {
            for j in i+1..n {
                get_mfe(&mut min_energy, &mut coax_energy, 
                        segments[i].len() - 2, i, 
                        segments[j].len() - 2, j);
            }
        }

        let ilen = segments[n-1].len();
        let jlen = segments[0].len();

        if jlen >= 2 {
                min_energy.get(0, n - 1, 1, 1)
                    .min(min_energy.get(0, n - 1, 1, 0))
                    .min(min_energy.get(0, n - 1, 0, 1))
                    .min(min_energy.get(0, n - 1, 0, 0))
        } else if jlen == 1 && ilen >= 2 {
                assert!(min_energy.get(0, n - 1, 1, 1) > i32::MAX/4); 
                min_energy.get(0, n - 1, 1, 1)
                    .min(min_energy.get(0, n - 1, 1, 0))
                    .min(min_energy.get(0, n - 1, 0, 1))
                    .min(min_energy.get(0, n - 1, 0, 0))
                    .min(min_energy.get(1, n - 2, 1, 1) + coax_energy.get(n-1, 0, 1, 1))
                    .min(min_energy.get(1, n - 2, 1, 1) + coax_energy.get(n-1, 0, 0, 1))
                    .min(min_energy.get(1, n - 2, 1, 0) + coax_energy.get(n-1, 0, 1, 1))
                    .min(min_energy.get(1, n - 2, 1, 0) + coax_energy.get(n-1, 0, 0, 1));
        }
        //println!("done");

        CoaxialStackingDP {
            segments,
            branching_pairs,
            min_energy,
            coax_energy,
            model,
        }
    }
    
    /// Forward pass filling E^{ab}[i..j]
    pub fn compute(&mut self) -> i32 {
        let n = self.branching_pairs.len();
                println!("{:?}", self.min_energy);
                println!("{:?}", self.coax_energy);
        for j in 1..n {
            for d5 in 0..=1 {
                for d3 in 0..=1 {
                    println!("i={} j={} d5={} d3={}", 0, j, d5, d3);
                    let val = self.update_e(0, j, d5, d3);
                    println!("val={}", val);
                    self.min_energy.set(0, j, d5, d3, val);
                }
            }
        }
        //for len in 2..=n {
        //    for i in 0..=n - len {
        //        let j = i + len - 1;
        //        for d5 in 0..=1 {
        //            for d3 in 0..=1 {
        //                println!("i={} j={} d5={} d3={}", i, j, d5, d3);
        //                let val = self.update_e(i, j, d5, d3);
        //                println!("val={}", val);
        //                self.min_energy.set(i, j, d5, d3, val);
        //            }
        //        }
        //        println!("{:?}", self.min_energy);
        //        println!("{:?}", self.coax_energy);
        //    }
        //}
        self.result()
    }

    fn update_e(&self, i: usize, j: usize, d5: usize, d3: usize) -> i32 {
        let mut best = i32::MAX / 2;

        // ----------------------------
        // Dangling contributions
        // ----------------------------
        let unpaired = self.segments[j].len() - 2;

        if unpaired >= 2 {
            best = best
                .min(self.min_energy.get(i, j - 1, d5, 1) + self.min_energy.get(j, j, 1, d3))
                .min(self.min_energy.get(i, j - 1, d5, 1) + self.min_energy.get(j, j, 0, d3))
                .min(self.min_energy.get(i, j - 1, d5, 0) + self.min_energy.get(j, j, 1, d3))
                .min(self.min_energy.get(i, j - 1, d5, 0) + self.min_energy.get(j, j, 0, d3));
        } else if unpaired == 1 {
            best = best
                .min(self.min_energy.get(i, j - 1, d5, 1) + self.min_energy.get(j, j, 0, d3))
                .min(self.min_energy.get(i, j - 1, d5, 0) + self.min_energy.get(j, j, 1, d3))
                .min(self.min_energy.get(i, j - 1, d5, 0) + self.min_energy.get(j, j, 0, d3));
            if j == 1 {
                best = best.min(self.coax_energy.get(j-1, j, d5, d3));
            } else {
                best = best.min(self.min_energy.get(i, j - 2, d5, 0) + self.coax_energy.get(j-1, j, 0, d3));
            }
        } else if unpaired == 0 {
            best = best
                .min(self.min_energy.get(i, j - 1, d5, 0) + self.min_energy.get(j, j, 0, d3));
            if j == 1 {
                best = best.min(self.coax_energy.get(j-1, j, d5, d3));
            } else {
                best = best.min(self.min_energy.get(i, j - 2, d5, 0) + self.coax_energy.get(j-1, j, 0, d3));
            }
        }
        best
    }

    pub fn from_exterior_loop(segments: &'a [&'a [Base]], model: &'a ViennaRNA) -> Self {
        let n = segments.len() - 1;
        let etables = &model.energy_tables;
        let mut branching_pairs = Vec::with_capacity(n);
        let mut min_energy = DPTable::new(n);
        let mut coax_energy = DPTable::new(n);
        for i in 0..n {
            let j = i + 1; 
            let pair = pair_type(*segments[i].last().unwrap(), segments[j][0]);
            println!("{i} {j} {:?}", pair);
            branching_pairs.push(pair);
            let d5 = d5_base(segments[i], i == 0);
            let d3 = d3_base(segments[j], j == n - 1);
            set_min_energy(&mut min_energy, i, d5, pair, d3, etables);

            if j >= n - 1 {
                continue;
            }

            let slen = segments[j].len() - 2;  
            if slen <= 1 {
                let j2 = i + 2; 
                let mm = if slen == 1 { segments[i].get(1) } else { None };
                let npair = pair_type(*segments[j].last().unwrap(), segments[j2][0]);
                set_coax_energy(&mut coax_energy, i, d5, pair, mm, npair, d3, model);
            }
        }

        CoaxialStackingDP {
            segments,
            branching_pairs,
            min_energy,
            coax_energy,
            model,
        }
    }
    
    /// Extract global minimum energy
    fn result(&self) -> i32 {
        let n = self.branching_pairs.len();
        let mut best = self.min_energy.get(0, n-1, 1, 1);
        for d5 in 0..=1 {
            for d3 in 0..=1 {
                best = best
                    .min(self.min_energy.get(0, n-1, d5, d3));
            }
        }
        best
    }
}

fn d5_base(seg: &[Base], exterior: bool) -> Option<&Base> {
    if exterior && seg.len() >= 2 {
        seg.get(seg.len() - 2)
    } else if !exterior && seg.len() >= 3 {
        seg.get(seg.len() - 2)
    } else {
        None
    }
}

fn d3_base(seg: &[Base], exterior: bool) -> Option<&Base> {
    if exterior && seg.len() >= 2 {
        seg.get(1)
    } else if !exterior && seg.len() >= 3 {
        seg.get(1)
    } else {
        None
    }
}

