/// acfp_domain_design.rs
///
/// Contains the struct Segment to map each position in the ACFP. 
/// Contains the DomainSequence (a wrapper for DomainRefVec).
///
///

use std::collections::VecDeque;
use rustc_hash::FxHashMap;
use nohash_hasher::{IntSet, IntMap};

use structure::{DotBracket, DotBracketVec, PairTable, StructureError};

use crate::domain::{DomainRef, DomainRefVec, DomainRegistry};
use crate::acfps::Acfp;
use crate::checks::dlfolding::{build_pair_scores, nussinov, traceback_structures};

#[derive(Debug)]
struct SegmentSequence {
    segments: Vec<Segment>,
}

impl SegmentSequence {
    /// To translate the DomainSequence from an ACFP, we also use a DomainRegistry,
    /// that can be used to specify Domain-lengths upfront. 
    pub fn design_from_acfp(acfp: &Acfp, registry: &mut DomainRegistry) -> Result<Self, String> {
        let pairh = acfp.pair_hierarchy().ok_or("Missing pair hierarchy")?;
        let ccomp = acfp.connected_components();
        let segments = design_segments(&ccomp, &pairh, registry)
            .into_iter()
            .map(|s| s.unwrap())
            .collect();

        Ok(SegmentSequence { segments })
    }
    
    pub fn get_detailed_segment_complexes(&self, registry: &DomainRegistry) -> () {
        let sequence = &self.get_domain_sequence();
        let p = build_pair_scores(sequence, &registry);
        let dp = nussinov(&p);
        
        for l in 0..sequence.len() {
            let all_structs = traceback_structures(0, l, &dp, &p);
            println!("[{}:] {}", l+1, sequence[0..=l].iter().map(|d| format!("{}", d)).collect::<Vec<_>>().join(" "));
            for st in all_structs {
                println!("{}", DotBracketVec::try_from(&PairTable(st)).unwrap());
            }
        }
    }

    pub fn get_completed_segment_complexes(&self, registry: &DomainRegistry) -> () {
        let sequence = &self.get_domain_sequence();
        let p = build_pair_scores(sequence, &registry);
        let dp = nussinov(&p);

        let mut lt = 0;
        for s in self.segments.iter() {
            lt += s.len();
            let all_structs = traceback_structures(0, lt-1, &dp, &p);
            println!("[{}:] {}", lt, sequence[0..lt].iter().map(|d| format!("{}", d)).collect::<Vec<_>>().join(" "));
            for st in all_structs {
                println!("{}", DotBracketVec::try_from(&PairTable(st)).unwrap());
            }
        }
    }

    pub fn get_segment_acfp(&self, registry: &DomainRegistry) -> Result<(), String> {
        let sequence = &self.get_domain_sequence();
        let p = build_pair_scores(sequence, registry);
        let dp = nussinov(&p);

        // Collect indices of all logic domains
        let logic_indices: std::collections::HashSet<usize> = self.segments.iter()
            .enumerate()
            .flat_map(|(si, seg)| {
                let offset: usize = self.segments[..si].iter().map(|s| s.len()).sum();
                std::iter::once(offset + seg.left_support.len()) // index of logic domain in full sequence
            })
        .collect();

        let mut lt = 0;
        for s in self.segments.iter() {
            lt += s.len();
            let all_structs = traceback_structures(0, lt - 1, &dp, &p);

            println!(
                "[{}:] {}",
                lt,
                sequence[0..lt]
                .iter()
                .map(|d| format!("{}", d))
                .collect::<Vec<_>>()
                .join(" ")
            );

            for st in all_structs {
                let mut dbv = vec![];
                
                for (i, &j_opt) in st.iter().enumerate() {
                    if !logic_indices.contains(&i) {
                        continue;
                    }
                    match j_opt {
                        Some(j) => {
                            if j > i { 
                                dbv.push(DotBracket::Open);
                            } else if j < i {
                                dbv.push(DotBracket::Close);
                            } else {
                                return Err(StructureError::InvalidPairTable(i).to_string());
                            }
                        },
                        None => dbv.push(DotBracket::Unpaired),
                    }
                }

                println!("{:?}", dbv);
            }
        }
        Ok(())
    }


    pub fn get_domain_sequence(&self) -> DomainRefVec {
        self.segments
            .iter()
            .flat_map(|seg| seg.full_sequence())
            .collect()
    }
        
}

#[derive(Clone, Debug)]
pub struct Segment {
    left_support: DomainRefVec,  // ["s1", "s2"], or ["s2*", "s1*"]
    logic: DomainRef,               // "L1", "L2*", etc.
    right_support: DomainRefVec, // ["s1", "s2"], or ["s2*", "s1*"]
    timer: DomainRef,
}

impl Segment {
    pub fn full_sequence(&self) -> DomainRefVec {
        self.left_support
            .iter()
            .cloned()
            .chain(std::iter::once(self.logic.clone()))
            .chain(self.right_support.iter().cloned())
            .chain(std::iter::once(self.timer.clone()))
            .collect()
    }

    pub fn left_support(&self) -> &DomainRefVec {
        &self.left_support
    }
    
    pub fn len(&self) -> usize {
        self.full_sequence().len()
    }
}

struct DomainBookkeeper {
    next_logic_id: usize,
    next_support_id: usize,
    support_length: usize,
    logic_length: usize,
}

impl DomainBookkeeper {
    pub fn new(support_length: usize, logic_length: usize) -> Self {
        Self {
            next_logic_id: 0,
            next_support_id: 0,
            support_length,
            logic_length,
        }
    }

    pub fn default() -> Self {
        Self {
            next_logic_id: 0,
            next_support_id: 0,
            support_length: 3,
            logic_length: 8,
        }
    }
 
    pub fn new_logic_pair(&mut self, i: usize, j: usize, weight: usize, registry: &mut DomainRegistry) -> (Segment, Segment) {
        self.next_logic_id += 1;

        let mut left_support = Vec::new();
        let mut right_support = Vec::new();

        for k in 0..(2 * weight) {
            self.next_support_id += 1;
            let name = format!("s{}", self.next_support_id);
            if k % 2 == 0 {
                left_support.insert(0, registry.intern(&name, 3));
            } else {
                right_support.push(registry.intern(&name, 3));
            }
        }

        (
            Segment {
                left_support: left_support.clone(),
                logic: registry.intern(&format!("L{}", self.next_logic_id), self.logic_length),
                right_support: right_support.clone(),
                timer: registry.intern(&format!("T{}", i), i),
            },
            Segment {
                left_support: right_support.into_iter().rev()
                    .map(|d| registry.get_complement(&d))
                    .collect(),
                logic: registry.intern(&format!("L{}*", self.next_logic_id), self.logic_length),
                right_support: left_support.into_iter().rev()
                    .map(|d| registry.get_complement(&d))
                    .collect(),
                timer: registry.intern(&format!("T{}", j), j),
            },
        )
    }

    fn extend_supports(&mut self, seg: &mut Segment, target_w: usize, registry: &mut DomainRegistry) {
        let comp = if seg.logic.name.ends_with('*') { "*".to_string() } else { "".to_string() };

        let needed = 2 * target_w;
        while seg.left_support.len() + seg.right_support.len() < needed {
            self.next_support_id += 1;
            let new = format!("s{}{}", self.next_support_id, comp);

            if seg.logic.name.ends_with('*') {
                if seg.left_support.len() <= seg.right_support.len() {
                    seg.left_support.insert(0, registry.intern(&new, self.support_length));
                } else {
                    seg.right_support.push(registry.intern(&new, self.support_length));
                }
            } else {
                if seg.left_support.len() > seg.right_support.len() {
                    seg.right_support.push(registry.intern(&new, self.support_length));
                } else {
                    seg.left_support.insert(0, registry.intern(&new, self.support_length));
                }
            }
        }
    }

}

fn design_segments(
    components: &[Vec<usize>],
    pair_hierarchy: &FxHashMap<(usize, usize), usize>,
    registry: &mut DomainRegistry, 
) -> Vec<Option<Segment>> {
    let max_index = components.iter().flatten().copied().max().unwrap_or(0);
    let mut seq: Vec<Option<Segment>> = vec![None; max_index + 1];
    let mut bookkeeper = DomainBookkeeper::default();

    // Build index: base -> component ID
    let mut base_to_component: IntMap<usize, usize> = IntMap::default();
    for (cid, comp) in components.iter().enumerate() {
        for &idx in comp {
            base_to_component.insert(idx, cid);
        }
    }

    let mut seen_components = IntSet::default();

    // Turn pairs into a queue for deferral
    let mut queue: VecDeque<((usize, usize), usize)> = pair_hierarchy
        .iter()
        .map(|(&(i, j), &w)| ((i, j), w))
        .collect();

    while let Some(((i, j), w)) = queue.pop_front() {
        let comp_id = base_to_component[&i];
        assert_eq!(base_to_component[&j], comp_id);

        let si = &seq[i];
        let sj = &seq[j];

        match (si, sj) {
            (None, None) => {
                if !seen_components.contains(&comp_id) {
                    let (seg_i, seg_j) = bookkeeper.new_logic_pair(i, j, w, registry);
                    seq[i] = Some(seg_i);
                    seq[j] = Some(seg_j);
                    seen_components.insert(comp_id);
                } else {
                    // Defer this pair until one side is known
                    queue.push_back(((i, j), w));
                }
            }
            (Some(si), None) => {
                let mut seg_i = si.clone();
                bookkeeper.extend_supports(&mut seg_i, w, registry);
                let seg_j = Segment {
                    left_support: seg_i.right_support.iter().take(w).rev()
                        .map(|d| registry.get_complement(&d))
                        .collect(),
                    logic: registry.get_complement(&seg_i.logic),
                    right_support: seg_i.left_support.iter().rev().take(w)
                        .map(|d| registry.get_complement(&d))
                        .collect(),
                    timer: registry.intern(&format!("T{}", j), j),
                };
                seq[i] = Some(seg_i);
                seq[j] = Some(seg_j);
            }
            (None, Some(sj)) => {
                let mut seg_j = sj.clone();
                bookkeeper.extend_supports(&mut seg_j, w, registry);
                let seg_i = Segment {
                    left_support: seg_j.right_support.iter().take(w).rev()
                        .map(|d| registry.get_complement(&d))
                        .collect(),
                    logic: registry.get_complement(&seg_j.logic),
                    right_support: seg_j.left_support.iter().rev().take(w)
                        .map(|d| registry.get_complement(&d))
                        .collect(),
                    timer: registry.intern(&format!("T{}", i), i),
                };
                seq[i] = Some(seg_i);
                seq[j] = Some(seg_j);
            }
            (Some(si), Some(sj)) => {
                let mut seg_i = si.clone();
                let mut seg_j = sj.clone();
                assert!(registry.are_complements(&seg_i.logic, &seg_j.logic));
                if seg_i.left_support.len() < seg_j.left_support.len() {
                    assert!(seg_i.right_support.len() < seg_j.right_support.len());
                    // TODO: assert that existing support matches before we overwrite it.
                    bookkeeper.extend_supports(&mut seg_j, w, registry);
                    if seg_i.left_support.len() + seg_i.right_support.len() < w {
                        seg_i = Segment {
                            left_support: seg_j.right_support.iter().take(w).rev()
                                .map(|d| registry.get_complement(&d))
                                .collect(),
                            logic: registry.get_complement(&seg_j.logic),
                            right_support: seg_j.left_support.iter().rev().take(w)
                                .map(|d| registry.get_complement(&d))
                                .collect(),
                            timer: registry.intern(&format!("T{}", i), i),
                        };
                    }
                } else if seg_j.left_support.len() < seg_i.left_support.len() {
                    assert!(seg_j.right_support.len() < seg_i.right_support.len());
                    // TODO: assert that existing support matches before we overwrite it.
                    bookkeeper.extend_supports(&mut seg_i, w, registry);
                    if seg_j.left_support.len() + seg_j.right_support.len() < w {
                        seg_j = Segment {
                            left_support: seg_i.right_support.iter().take(w).rev()
                                .map(|d| registry.get_complement(&d))
                                .collect(),
                            logic: registry.get_complement(&seg_i.logic),
                            right_support: seg_i.left_support.iter().rev().take(w)
                                .map(|d| registry.get_complement(&d))
                                .collect(),
                            timer: registry.intern(&format!("T{}", j), j),
                        };
                    }
                }
                seq[i] = Some(seg_i);
                seq[j] = Some(seg_j);
            }
        }
    }

    seq[1..].to_vec()
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_design_segments_balancing_and_expansion() {
        // Define the pair hierarchy
        let pair_hierarchy: FxHashMap<(usize, usize), usize> = vec![
            ((1, 2), 1),
            ((3, 4), 2),
            ((2, 3), 3),
        ].into_iter().collect();

        // One connected component
        let components = vec![vec![1, 2, 3, 4]];

        // Call the function
        let mut registry = DomainRegistry::new();

        let segments = design_segments(&components, &pair_hierarchy, &mut registry);
        let sequence: DomainRefVec = segments.iter().flatten().flat_map(|d| d.full_sequence()).collect();
        println!("{:?}", sequence
                .iter()
                .map(|d| d.name.clone()).collect::<Vec<_>>()
                .join(" ")
            );

        // Print results for debugging
        for seg in segments.iter().flatten() {
            println!("{:?}", seg.full_sequence()
                .iter()
                .map(|d| d.name.clone()).collect::<Vec<_>>()
                .join(" "));
        }

        // Check symmetry and logic domain names
        let s1 = segments[0].as_ref().unwrap();
        let s2 = segments[1].as_ref().unwrap();
        let s3 = segments[2].as_ref().unwrap();
        let s4 = segments[3].as_ref().unwrap();

        // Ensure logic domains are complementary
        assert_eq!(s1.logic.name, "L1");
        assert_eq!(s2.logic.name, "L1*");
        assert_eq!(s3.logic.name, "L1");
        assert_eq!(s4.logic.name, "L1*");

        // Ensure left + right support lengths = 2w
        assert_eq!(s1.left_support.len() + s1.right_support.len(), 2);
        assert_eq!(s2.left_support.len() + s2.right_support.len(), 6);
        assert_eq!(s3.left_support.len() + s3.right_support.len(), 6);
        assert_eq!(s4.left_support.len() + s4.right_support.len(), 4);
    }

    #[test]
    fn test_design_acfp_04() {
        let acfp = Acfp::try_from(". () (). ()() (().) ").unwrap();
        assert!(acfp.is_valid());

        let mut registry = DomainRegistry::new();
        let segments = SegmentSequence::design_from_acfp(&acfp, &mut registry).unwrap();
        println!("{:?}", segments);

        let sequence = &segments.get_domain_sequence();
        println!("{}", sequence.iter().map(|d| format!("{}", d)).collect::<Vec<_>>().join(" "));

        let p = build_pair_scores(sequence, &registry);
        let dp = nussinov(&p);
        
        let mut lt = 0;
        for s in segments.segments.iter() {
            lt += s.len();
            let all_structs = traceback_structures(0, lt-1, &dp, &p);
            println!("{:?} {:?}", s.full_sequence().iter().map(|d| format!("{}", d)).collect::<Vec<_>>().join(" "), s.len()-1);
            println!("{}", DotBracketVec::try_from(&PairTable(all_structs[0].clone())).unwrap());
            assert_eq!(all_structs.len(), 1);
        }

        let _ = &segments.get_detailed_segment_complexes(&registry);
        println!();
        let _ = &segments.get_completed_segment_complexes(&registry);
        println!();
        let _ = &segments.get_segment_acfp(&registry);
        assert!(false);
    }


}
