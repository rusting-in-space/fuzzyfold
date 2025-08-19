
use structure::PairTable;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum NearestNeighborLoop {
    Hairpin {
        closing: (usize, usize), // (i, j)
    },
    Interior {
        closing: (usize, usize),
        inner: (usize, usize),
    },
    Multibranch {
        closing: (usize, usize),
        //NOTE: this list must ALWAYS be in 5'->3' order.
        branches: Vec<(usize, usize)>,
    },
    Exterior {
        //NOTE: this list must ALWAYS be in 5'->3' order.
        branches: Vec<(usize, usize)>,
    },
}

impl NearestNeighborLoop {
    pub fn classify(
        closing: Option<(usize, usize)>, 
        branches: Vec<(usize, usize)>, 
    ) -> Self {
        match closing {
            None => Self::Exterior { branches },
            Some((i, j)) => match branches.len() {
                0 => Self::Hairpin { closing: (i, j) },
                1 => Self::Interior { closing: (i, j), inner: branches[0] },
                _ => Self::Multibranch { closing: (i, j), branches },
            },
        }
    }
}

pub trait LoopDecomposition {
    fn for_each_loop<F: FnMut(&NearestNeighborLoop)>(&self, f: F);

    fn loops(&self) -> Vec<NearestNeighborLoop> {
        let mut out = Vec::new();
        self.for_each_loop(|l| out.push(l.clone()));
        out
    }

    fn loop_enclosed_by(&self, closing: Option<(usize, usize)>) -> NearestNeighborLoop;

    fn get_enclosing_pair(&self, i: usize, j: usize) -> Option<(usize, usize)>;
}

impl LoopDecomposition for PairTable {
    fn for_each_loop<F: FnMut(&NearestNeighborLoop)>(&self, mut f: F) {
        fn recurse<F: FnMut(&NearestNeighborLoop)>(
            pt: &PairTable,
            closing: Option<(usize, usize)>,
            f: &mut F,
        ) {
            let mut branches = Vec::new();

            let (mut p, j) = if let Some((i, j)) = closing {
                (i + 1, j) 
            } else { 
                (0, pt.len())
            };

            while p < j {
                if let Some(q) = pt[p] {
                    assert!(q > p);
                    branches.push((p, q));
                    // Recurse into child loop
                    recurse(pt, Some((p, q)), f);
                    p = q + 1;
                } else {
                    p += 1;
                }
            }
            f(&NearestNeighborLoop::classify(closing, branches));
        }
        recurse(&self, None, &mut f);
    }

    fn loop_enclosed_by(&self, closing: Option<(usize, usize)>) -> NearestNeighborLoop {
        let mut branches = Vec::new();

        let (mut p, j) = if let Some((i, j)) = closing {
            (i + 1, j) 
        } else { 
            (0, self.len())
        };

        while p < j {
            if let Some(q) = self[p] {
                assert!(q > p);
                branches.push((p, q));
                p = q + 1;
            } else {
                p += 1;
            }
        }
        NearestNeighborLoop::classify(closing, branches)
    }

    fn get_enclosing_pair(&self, i: usize, j: usize) -> Option<(usize, usize)> { 
        for q in j..self.len() {
            if let Some(p) = self[q] {
                if p < i {
                    return Some((p, q));
                }
            }
        }
        return None
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decompose_loops_empty() {
        let dbn = "......."; // all unpaired → exterior loop only
        let eloop = NearestNeighborLoop::Exterior { 
            branches: vec![], 
        };

        let loops = PairTable::try_from(dbn).expect("valid").loops();
        println!("{:?}", loops);
        assert_eq!(loops, vec![eloop]);
    }

    #[test]
    fn test_decompose_loops_hairpin() {
        let dbn = ".(...).";

        let eloop = NearestNeighborLoop::Exterior { 
            branches: vec![(1, 5)], 
        };

        let hloop = NearestNeighborLoop::Hairpin { 
            closing: (1, 5) 
        };

        let loops = PairTable::try_from(dbn).expect("valid").loops();
        println!("{:?}", loops);
        assert!(loops.len() == 2);
        assert!(loops.contains(&eloop));
        assert!(loops.contains(&hloop));
    }

    #[test]
    fn test_decompose_loops_wild() {
        let dbn = ".(.((...))()((()))).((...()))";
        let loops = PairTable::try_from(dbn).expect("valid").loops();
        println!("{:?}", loops);
    }
}

