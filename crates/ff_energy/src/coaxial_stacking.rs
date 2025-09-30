
use std::fmt;

use crate::Base;
use crate::PairTypeRNA;
use crate::EnergyTables;

const ALL_DANGLE_CONTEXTS: &[(usize, usize)] = &[(0, 0), (0, 1), (1, 0), (1, 1)];
const EXC_DANGLE_CONTEXTS: &[(usize, usize)] = &[(0, 0), (0, 1), (1, 0)];
const NON_DANGLE_CONTEXTS: &[(usize, usize)] = &[(0, 0)];

/// Unified structure for min_energy and coaxial_energy DP tables.
pub struct CoaxialStackingDP {
    n: usize,
    min_energy: Vec<i32>,
    coax_energy: Vec<i32>,
}

impl fmt::Debug for CoaxialStackingDP {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "CoaxialStackingDP (n = {})", self.n)?;
        for i in 0..self.n {
            for j in 0..self.n {
                let vals: Vec<String> = (0..=1)
                    .flat_map(|d5| (0..=1).map(move |d3| (d5, d3)))
                    .map(|(d5, d3)| {
                        let v = self.get_min(i, j, d5, d3);
                        if v >= i32::MAX / 4 {
                            format!("({},{}): ∞", d5, d3)
                        } else {
                            format!("({},{}): {}", d5, d3, v)
                        }
                    })
                    .collect();
                writeln!(f, "[{},{}] min -> {}", i, j, vals.join(", "))?;

                let vals: Vec<String> = (0..=1)
                    .flat_map(|d5| (0..=1).map(move |d3| (d5, d3)))
                    .map(|(d5, d3)| {
                        let v = self.get_coax(i, j, d5, d3);
                        if v >= i32::MAX / 4 {
                            format!("({},{}): ∞", d5, d3)
                        } else {
                            format!("({},{}): {}", d5, d3, v)
                        }
                    })
                    .collect();
                writeln!(f, "[{},{}] coax -> {}", i, j, vals.join(", "))?;
            }
        }
        Ok(())
    }
}

impl CoaxialStackingDP {
    pub fn new(n: usize) -> Self {
        Self {
            n,
            min_energy: vec![i32::MAX / 2; n * n * 4],
            coax_energy: vec![i32::MAX / 2; n * n * 4],
        }
    }

    fn idx(&self, i: usize, j: usize, d5: usize, d3: usize) -> usize {
        (((i * self.n) + j) << 2) + ((d5 << 1) | d3)
    }

    // ---- min energy accessors
    fn get_min(&self, i: usize, j: usize, d5: usize, d3: usize) -> i32 {
        self.min_energy[self.idx(i, j, d5, d3)]
    }
    fn set_min(&mut self, i: usize, j: usize, d5: usize, d3: usize, val: i32) {
        let idx = self.idx(i, j, d5, d3);
        self.min_energy[idx] = val;
    }

    // ---- coax energy accessors
    fn get_coax(&self, i: usize, j: usize, d5: usize, d3: usize) -> i32 {
        self.coax_energy[self.idx(i, j, d5, d3)]
    }
    fn set_coax(&mut self, i: usize, j: usize, d5: usize, d3: usize, val: i32) {
        let idx = self.idx(i, j, d5, d3);
        self.coax_energy[idx] = val;
    }

    pub fn compute_exterior_mfe(&mut self, segments: &[&[Base]]) -> i32 {
        assert_eq!(self.n, segments.len() - 1);
        for j in 1..self.n {
            self.get_mfe(
                segments[0].len() - 2, 0,
                segments[j].len() - 2, j,
                self.n,
            );
        }

        [
            self.get_min(0, self.n-1, 1, 1),
            self.get_min(0, self.n-1, 1, 0),
            self.get_min(0, self.n-1, 0, 1),
            self.get_min(0, self.n-1, 0, 0),
        ]
        .into_iter()
        .min()
        .unwrap_or(i32::MAX)
    }

    pub fn compute_multibranch_mfe(&mut self, segments: &[&[Base]]) -> i32 {
        assert_eq!(self.n, segments.len());
        for i in 0..=1 {
            for j in 1..self.n {
                let j = (i + j) % self.n;
                self.get_mfe(
                    segments[i].len() - 2, i, 
                    segments[j].len() - 2, j,
                    self.n
                );
            }
        }
 
        [
            self.get_min(0, self.n-1, 1, 1),
            self.get_min(0, self.n-1, 1, 0),
            self.get_min(0, self.n-1, 0, 1),
            self.get_min(0, self.n-1, 0, 0),
            self.get_min(1, 0, 1, 1),
            self.get_min(1, 0, 1, 0),
            self.get_min(1, 0, 0, 1),
            self.get_min(1, 0, 0, 0),
        ]
        .into_iter()
        .min()
        .unwrap_or(i32::MAX)
    }

    // ---- your methods (logic unchanged, only adapted to unified struct)
    pub fn set_min_energy(&mut self,
        i: usize,
        d5: Option<&Base>,
        pi: PairTypeRNA,
        d3: Option<&Base>,
        etables: &EnergyTables,
    ) {
        self.set_min(i, i, 0, 0, 0);
        self.set_min(i, i, 1, 0,
            d5.map_or(i32::MAX/2, |&b| {
                etables.dangle5[pi as usize][b as usize].unwrap()
            })
        );
        self.set_min(i, i, 0, 1,
            d3.map_or(i32::MAX/2, |&b| {
                etables.dangle3[pi as usize][b as usize].unwrap()
            })
        );
        self.set_min(i, i, 1, 1,
            match (d5, d3) {
                (Some(&b5), Some(&b3)) => etables
                    .mismatch_exterior[pi as usize][b5 as usize][b3 as usize].unwrap(),
                _ => i32::MAX/2,
            }
        );
    }

    fn set_min_energy_combinations(&mut self,
        i: usize,
        j: usize,
        outer_context: &[(usize, usize)],
        inner_context: &[(usize, usize)],
    ) {
        let pj = wrap_sub1(j, self.n);
        for &(d5, d3) in outer_context {
            let mut best = i32::MAX;
            for &(id5, id3) in inner_context {
                let val = self.get_min(i, pj, d5, id3) + self.get_min(j, j, id5, d3);
                best = best.min(val);
            }
            self.set_min(i, j, d5, d3, best);
        }
    }

    pub fn set_coax_energy(&mut self,
        i: usize,
        j: usize,
        d5: Option<&Base>,
        pi: PairTypeRNA,
        mm: Option<&Base>,
        pj: PairTypeRNA,
        d3: Option<&Base>,
        etables: &EnergyTables,
        coaxial_mm_discontious_en37: i32,
    ) {
        if let Some(&m) = mm {
            self.set_coax(i, j, 1, 0,
                d5.map_or(i32::MAX/2, |&b| {
                    let pair = PairTypeRNA::from((b, m));
                    let bonus = if pair.is_wobble() {
                        -20
                    } else if pair.can_pair() {
                        -40
                    } else { 0 };
                    etables.mismatch_exterior[pi as usize][b as usize][m as usize].unwrap()
                        + coaxial_mm_discontious_en37 + bonus
                }));
            self.set_coax(i, j, 0, 1,
                d3.map_or(i32::MAX/2, |&b| {
                    let pair = PairTypeRNA::from((m, b));
                    let bonus = if pair.is_wobble() {
                        -20
                    } else if pair.can_pair() {
                        -40
                    } else { 0 };
                    etables.mismatch_exterior[pi as usize][m as usize][b as usize].unwrap()
                        + coaxial_mm_discontious_en37 + bonus
                }));
        } else {
            self.set_coax(i, j, 0, 0, {
                let pj = pj.invert();
                etables.stack[pi as usize][pj as usize].unwrap()
            });
        }
    }

    fn set_coax_energy_combinations(&mut self,
        i: usize,
        j: usize,
        outer_context: &[(usize, usize)],
        inner_context: &[(usize, usize)],
        coax_context: &[(usize, usize)],
    ) {
        let pj = wrap_sub1(j, self.n);
        let ppj = wrap_sub1(pj, self.n);
        for &(d5, d3) in outer_context {
            let mut best = i32::MAX;
            for &(id5, id3) in inner_context {
                let val = self.get_min(i, pj, d5, id3) + self.get_min(j, j, id5, d3);
                best = best.min(val);
            }
            for &(id5, id3) in coax_context {
                let val = self.get_min(i, ppj, d5, id3) + self.get_coax(pj, j, id5, d3);
                best = best.min(val);
            }
            self.set_min(i, j, d5, d3, best);
        }
    }

    fn set_coax_init_combinations(&mut self,
        i: usize,
        j: usize,
        outer_context: &[(usize, usize)],
        inner_context: &[(usize, usize)],
    ) {
        let pj = wrap_sub1(j, self.n);
        for &(d5, d3) in outer_context {
            let mut best = i32::MAX;
            for &(id5, id3) in inner_context {
                let val = self.get_min(i, pj, d5, id3) + self.get_min(j, j, id5, d3);
                best = best.min(val);
            }
            best = best.min(self.get_coax(i, j, d5, d3));
            self.set_min(i, j, d5, d3, best);
        }
    }

    fn get_mfe(&mut self,
        ilen: usize,
        i: usize,
        jlen: usize,
        j: usize,
        n: usize,
    ) {
        match jlen {
            0 => {
                if j == wrap_add1(i, n) {
                    self.set_coax_init_combinations(i, j,
                        ALL_DANGLE_CONTEXTS,
                        NON_DANGLE_CONTEXTS)
                } else if ilen >= 2 {
                    self.set_coax_energy_combinations(i, j,
                        ALL_DANGLE_CONTEXTS,
                        NON_DANGLE_CONTEXTS,
                        ALL_DANGLE_CONTEXTS)
                } else if ilen == 1 {
                    self.set_coax_energy_combinations(i, j,
                        ALL_DANGLE_CONTEXTS,
                        NON_DANGLE_CONTEXTS,
                        EXC_DANGLE_CONTEXTS)
                } else {
                    assert!(ilen == 0);
                    self.set_coax_energy_combinations(i, j,
                        ALL_DANGLE_CONTEXTS,
                        NON_DANGLE_CONTEXTS,
                        NON_DANGLE_CONTEXTS)
                }
            },
            1 => {
                if j == wrap_add1(i, n) {
                    self.set_coax_init_combinations(i, j,
                        ALL_DANGLE_CONTEXTS,
                        EXC_DANGLE_CONTEXTS)
                } else if ilen >= 2 {
                    self.set_coax_energy_combinations(i, j,
                        ALL_DANGLE_CONTEXTS,
                        EXC_DANGLE_CONTEXTS,
                        ALL_DANGLE_CONTEXTS)
                } else if ilen == 1 {
                    self.set_coax_energy_combinations(i, j,
                        ALL_DANGLE_CONTEXTS,
                        EXC_DANGLE_CONTEXTS,
                        EXC_DANGLE_CONTEXTS)
                } else {
                    assert!(ilen == 0);
                    self.set_coax_energy_combinations(i, j,
                        ALL_DANGLE_CONTEXTS,
                        EXC_DANGLE_CONTEXTS,
                        NON_DANGLE_CONTEXTS)
                }
            },
            _ => {
                self.set_min_energy_combinations(i, j,
                    ALL_DANGLE_CONTEXTS, ALL_DANGLE_CONTEXTS)
            }
        }
    }
}
 
fn wrap_sub1(i: usize, n: usize) -> usize {
    // Avoiding modulo, which may be slow.
    if i > 0 { i - 1 } else { n - 1 }
}

fn wrap_add1(i: usize, n: usize) -> usize {
    // Avoiding modulo, which may be slow.
    if i + 1 < n { i + 1 } else { 0 }
}

