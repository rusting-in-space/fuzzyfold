
use ff_structure::PairTable;
use crate::NearestNeighborLoop;

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
        recurse(self, None, &mut f);
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
        None
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
    fn test_loop_ranges() {
        let hloop = NearestNeighborLoop::Hairpin { 
            closing: (1, 5) 
        };
        println!("{:?}", hloop.unpaired_indices(26));
        println!("{:?}", hloop.inclusive_unpaired_indices(26));

        let iloop = NearestNeighborLoop::Interior { 
            closing: (1, 9),
            inner: (2, 8) 
        };
        println!("{:?}", iloop.unpaired_indices(26));
        println!("{:?}", iloop.inclusive_unpaired_indices(26));
 
        let iloop = NearestNeighborLoop::Interior { 
            closing: (1, 9),
            inner: (2, 7) 
        };
        println!("{:?}", iloop.unpaired_indices(26));
        println!("{:?}", iloop.inclusive_unpaired_indices(26));
 

        let mloop = NearestNeighborLoop::Multibranch { 
            closing: (1, 15),
            branches: vec![(2, 4), (5, 9)],
        };
        println!("{:?}", mloop.unpaired_indices(26));
        println!("{:?}", mloop.inclusive_unpaired_indices(26));

        let eloop = NearestNeighborLoop::Exterior { 
            branches: vec![(1, 5), (6, 11)], 
        };
        println!("{:?}", eloop.unpaired_indices(26));
        println!("{:?}", eloop.inclusive_unpaired_indices(26));

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

