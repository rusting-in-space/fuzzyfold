
use colored::*; // brings in the `.red()`, `.blue()` etc.
use std::fmt;
use std::ops::Range;
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

impl fmt::Display for NearestNeighborLoop {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NearestNeighborLoop::Hairpin { closing: (i, j) } => {
                write!(f, "{} ({}, {})", "Hairpin loop".red().bold(), i, j)
            }
            NearestNeighborLoop::Interior { closing: (i, j), inner: (p, q) } => {
                write!(f, "{} ({}, {}), ({}, {})", "Interior loop".blue().bold(), i, j, p, q)
            }
            NearestNeighborLoop::Multibranch { closing: (i, j), branches } => {
                write!(f, "{} ({}, {}), ({} branches)", "Multibranch".green().bold(), i, j, branches.len())
            }
            NearestNeighborLoop::Exterior { branches } => {
                write!(f, "{} ({} branches)", "Exterior loop".cyan().bold(), branches.len())
            }
        }
    }
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

    pub fn closing(&self) -> Option<(usize, usize)> {
        match self { Self::Hairpin { closing }
            | Self::Interior { closing, .. }
            | Self::Multibranch { closing, .. } => Some(*closing),
            Self::Exterior { .. } => None,
        }
    }

    fn unpaired_ranges(&self, len: usize) -> Vec<Range<usize>> {
        match self {
            Self::Hairpin { closing: (i, j) } => { 
                vec![(*i + 1..*j)]
            },
            Self::Interior { closing: (i, j),  inner: (p, q) } => { 
                vec![(*i + 1)..(*p), (*q + 1)..(*j)]
            },
            Self::Multibranch { closing: (i, j), branches } => {
                let mut result = vec![];
                let mut start = *i;
                for &(p, q) in branches {
                    result.push((start+1)..p);
                    start = q;
                }
                result.push((start+1)..(*j));
                result
            }
            Self::Exterior { branches } => {
                let mut result = Vec::new();
                let mut start = 0;
                for &(p, q) in branches {
                    result.push(start..p);
                    start = q+1;
                }
                result.push(start..len);
                result
            }
        }
    }

    pub fn unpaired_indices(&self, len: usize) -> Vec<usize> {
        self.unpaired_ranges(len)
            .into_iter()
            .flat_map(|r| r.collect::<Vec<_>>())
            .collect()
    }

    /// Returns all sequence indices that should point to this loop.
    pub fn inclusive_unpaired_indices(&self, len: usize) -> Vec<usize> {
        self.unpaired_ranges(len-1)
            .into_iter()
            .map(|r| r.start..=r.end)
            .flat_map(|r| r.collect::<Vec<_>>())
            .collect()
    }
   
    /// Split the given loop into two new loops at the indices i,j
    /// NOTE: Returns (outer, inner)
    pub fn split_loop(&self, i: usize, j: usize) -> (Self, Self) {
        assert!(i < j, "Split pair (i,j) must satisfy i < j");
        match self {
            Self::Hairpin { closing: (a, b) } => {
                assert!(*a < i && j < *b, "Pair (i,j) must be within hairpin loop");
                (Self::Interior { closing: (*a, *b), inner: (i, j), },
                 Self::Hairpin { closing: (i, j), })
            }
            
            Self::Interior { closing: (a, b), inner: (p, q) } => {
                assert!(*a < i && j < *b, "Pair (i,j) outside of loop");
                assert!(!(*p < i && j < *q), "Pair (i,j) outside of loop");

                if i < *p && *q < j {
                    (Self::Interior { closing: (*a, *b), inner: (i, j) },
                     Self::Interior { closing: (i, j), inner: (*p, *q) })
                } else if j < *p {
                    (Self::Multibranch { closing: (*a, *b), branches: vec![(i, j), (*p, *q)] },
                     Self::Hairpin { closing: (i, j) })
                } else if *q < i {
                    (Self::Multibranch { closing: (*a, *b), branches: vec![(*p, *q), (i, j)] },
                     Self::Hairpin { closing: (i, j) })
                } else {
                    panic!("that really should not happen.");
                }
            }

            Self::Multibranch { closing: (a, b), branches } => {
                assert!(*a < i && j < *b, "Pair (i,j) outside loop");

                let mut outer_branches = vec![(i, j)];
                let mut inner_branches = vec![];
    
                for &(p, q) in branches {
                    assert!(p < q);
                    if j < p || q < i {
                        outer_branches.push((p, q));
                    } else { 
                        assert!(i < p && q < j);
                        inner_branches.push((p, q));
                    } 
                }

                outer_branches.sort_unstable();
                inner_branches.sort_unstable();
                (Self::classify(Some((*a, *b)), outer_branches),
                 Self::classify(Some((i, j)), inner_branches))
            }

            Self::Exterior { branches } => {
                let mut outer_branches = vec![];
                let mut inner_branches = vec![];
    
                outer_branches.push((i, j));
                for &(p, q) in branches {
                    assert!(p < q);
                    if j < p || q < i {
                        outer_branches.push((p, q));
                    } else { 
                        assert!(i < p && q < j);
                        inner_branches.push((p, q));
                    } 
                }

                outer_branches.sort_unstable();
                inner_branches.sort_unstable();
                (Self::classify(None, outer_branches),
                 Self::classify(Some((i, j)), inner_branches))
            }
        }
    }

    /// Join two loops by reference. 
    /// NOTE: the outer loop must be joined with the inner loop!
    pub fn join_loop(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Hairpin { .. }, _) => { 
                panic!("A hairpin cannot be the outer loop!");
            }

            (_, Self::Exterior { .. }) => { 
                panic!("Multi-stranded moves are not supported.");
            }

            (Self::Interior { closing: outer_closing, inner },
             Self::Hairpin { closing: inner_closing }) => {
                assert_eq!(inner, inner_closing, "Cannot join interior & haipin loops!");
                Self::Hairpin { closing: *outer_closing } 
            }

            (Self::Interior { closing: outer_closing, inner: outer_inner },
             Self::Interior { closing: inner_closing, inner: inner_inner })  => {
                assert_eq!(outer_inner, inner_closing, "Cannot join interior & interior loops!");
                Self::Interior { closing: *outer_closing, inner: *inner_inner }
            }
 
            (Self::Interior { closing: outer_closing, inner },
             Self::Multibranch { closing: inner_closing, branches }) => {
                assert_eq!(inner, inner_closing, "Cannot join interior & multibranch loops!");
                Self::Multibranch { closing: *outer_closing, branches: branches.clone() }
            },

            (Self::Multibranch { closing: outer_closing, branches },
             Self::Hairpin { closing: inner_closing }
            ) => {
                let new_branches: Vec<_> = branches.iter().cloned()
                    .filter(|b| b != inner_closing)
                    .collect();
                assert_eq!(branches.len(), new_branches.len() + 1, "Cannot join multibranch & hairpin loops!");
                Self::classify(Some(*outer_closing), new_branches)
            }
  
            (Self::Multibranch { closing: outer_closing, branches },
             Self::Interior { closing: inner_closing, inner }
            ) => {
                let new_branches: Vec<_>  = branches.iter().cloned()
                    .map(|x| if x == *inner_closing { *inner } else { x })
                    .collect();
                Self::Multibranch { closing: *outer_closing, branches: new_branches }
            },
 
            (Self::Multibranch { closing: outer_closing, branches: outer_branches },
             Self::Multibranch { closing: inner_closing, branches: inner_branches }
            ) => {
                let new_branches: Vec<_>  = outer_branches.iter().cloned()
                    .flat_map(|x| if x == *inner_closing { inner_branches.clone()  } else { vec![x] })
                    .collect();
                Self::Multibranch { closing: *outer_closing, branches: new_branches }
            },
            
            (Self::Exterior { branches },
             Self::Hairpin { closing: inner_closing }
            ) => {
                let new_branches: Vec<_>  = branches.iter().cloned()
                    .filter(|b| b != inner_closing)
                    .collect();
                Self::Exterior { branches: new_branches }
            }

            (Self::Exterior { branches },
             Self::Interior { closing: inner_closing, inner }
            ) => {
                let new_branches: Vec<_>  = branches.iter().cloned()
                    .map(|x| if x == *inner_closing { *inner } else { x })
                    .collect();
                Self::Exterior { branches: new_branches }
            },
              
            (Self::Exterior { branches: outer_branches },
             Self::Multibranch { closing: inner_closing, branches: inner_branches }
            ) => {
                let new_branches: Vec<_>  = outer_branches.iter().cloned()
                    .flat_map(|x| if x == *inner_closing { inner_branches.clone()  } else { vec![x] })
                    .collect();
                Self::Exterior { branches: new_branches }
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
    fn test_loop_ranges() {
        let hloop = NearestNeighborLoop::Hairpin { 
            closing: (1, 5) 
        };
        println!("{:?}", hloop.unpaired_ranges(26));
        println!("{:?}", hloop.unpaired_indices(26));
        println!("{:?}", hloop.inclusive_unpaired_indices(26));

        let iloop = NearestNeighborLoop::Interior { 
            closing: (1, 9),
            inner: (2, 8) 
        };
        println!("{:?}", iloop.unpaired_ranges(26));
        println!("{:?}", iloop.unpaired_indices(26));
        println!("{:?}", iloop.inclusive_unpaired_indices(26));
 
        let iloop = NearestNeighborLoop::Interior { 
            closing: (1, 9),
            inner: (2, 7) 
        };
        println!("{:?}", iloop.unpaired_ranges(26));
        println!("{:?}", iloop.unpaired_indices(26));
        println!("{:?}", iloop.inclusive_unpaired_indices(26));
 

        let mloop = NearestNeighborLoop::Multibranch { 
            closing: (1, 15),
            branches: vec![(2, 4), (5, 9)],
        };
        println!("{:?}", mloop.unpaired_ranges(26));
        println!("{:?}", mloop.unpaired_indices(26));
        println!("{:?}", mloop.inclusive_unpaired_indices(26));

        let eloop = NearestNeighborLoop::Exterior { 
            branches: vec![(1, 5), (6, 11)], 
        };
        println!("{:?}", eloop.unpaired_ranges(26));
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

