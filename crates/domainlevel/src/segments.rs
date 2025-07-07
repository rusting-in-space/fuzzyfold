
use std::collections::VecDeque;
use rustc_hash::FxHashMap;
use nohash_hasher::{IntSet, IntMap};

#[derive(Clone, Debug)]
pub struct Segment {
    left_support: Vec<String>,  // ["s1", "s2"], or ["s2*", "s1*"]
    logic: String,             // "L1", "L2*", etc.
    right_support: Vec<String>, // ["s1", "s2"], or ["s2*", "s1*"]
    timer: String,
}

impl Segment {
    pub fn full_sequence(&self) -> Vec<String> {
        self.left_support
            .iter()
            .cloned()
            .chain(std::iter::once(self.logic.clone()))
            .chain(self.right_support.iter().cloned())
            .chain(std::iter::once(self.timer.clone()))
            .collect()
    }

    pub fn left_support(&self) -> &Vec<String> {
        &self.left_support
    }
}

#[derive(Default)]
struct DomainBookkeeper {
    next_logic_id: usize,
    next_support_id: usize,
}

impl DomainBookkeeper {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn new_logic_pair(&mut self, i: usize, j: usize, weight: usize) -> (Segment, Segment) {
        self.next_logic_id += 1;

        let mut left_support = Vec::new();
        let mut right_support = Vec::new();

        for k in 0..(2 * weight) {
            self.next_support_id += 1;
            let name = format!("s{}", self.next_support_id);
            if k % 2 == 0 {
                left_support.insert(0, name);
            } else {
                right_support.push(name);
            }
        }

        (
            Segment {
                left_support: left_support.clone(),
                logic: format!("L{}", self.next_logic_id),
                right_support: right_support.clone(),
                timer: format!("T{}", i),
            },
            Segment {
                left_support: right_support.into_iter().rev().map(|s| format!("{}*", s)).collect(),
                logic: format!("L{}*", self.next_logic_id),
                right_support: left_support.into_iter().rev().map(|s| format!("{}*", s)).collect(),
                timer: format!("T{}", j),
            },
        )
    }

}

fn extend_supports(
    seg: &mut Segment,
    bk: &mut DomainBookkeeper,
    target_w: usize,
) {
    let comp = if seg.logic.ends_with('*') { "*".to_string() } else { "".to_string() };
    let needed = 2 * target_w;
    while seg.left_support.len() + seg.right_support.len() < needed {
        bk.next_support_id += 1;
        let new = format!("s{}{}", bk.next_support_id, comp);
        if seg.logic.ends_with('*') {
            if seg.left_support.len() <= seg.right_support.len() {
                seg.left_support.insert(0, new);
            } else {
                seg.right_support.push(new);
            }
        } else {
            if seg.left_support.len() > seg.right_support.len() {
                seg.right_support.push(new);
            } else {
                seg.left_support.insert(0, new);
            }
        }
    }
}


pub fn assign_segments(
    components: &[Vec<usize>],
    pair_hierarchy: &FxHashMap<(usize, usize), usize>,
) -> Vec<Option<Segment>> {
    let max_index = components.iter().flatten().copied().max().unwrap_or(0);
    let mut seq: Vec<Option<Segment>> = vec![None; max_index + 1];
    let mut bookkeeper = DomainBookkeeper::new();

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
        println!("Processing: {} {} {}", i, j, w);
        let comp_id = base_to_component[&i];
        assert_eq!(base_to_component[&j], comp_id);

        let si = &seq[i];
        let sj = &seq[j];

        match (si, sj) {
            (None, None) => {
                if !seen_components.contains(&comp_id) {
                    let (seg_i, seg_j) = bookkeeper.new_logic_pair(i, j, w);
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
                extend_supports(&mut seg_i, &mut bookkeeper, w);

                let comp = if seg_i.logic.ends_with('*') { "".to_string() } else { "*".to_string() };
                let seg_j = Segment {
                    left_support: seg_i.right_support.iter().take(w).rev().map(|s| format!("{}{}", s.trim_end_matches('*'), comp)).collect(),
                    logic: format!("{}{}", seg_i.logic.trim_end_matches('*'), comp),
                    right_support: seg_i.left_support.iter().rev().take(w).map(|s| format!("{}{}", s.trim_end_matches('*'), comp)).collect(),
                    timer: format!("T{}", j),
                };
                seq[i] = Some(seg_i);
                seq[j] = Some(seg_j);
            }
            (None, Some(sj)) => {
                let mut seg_j = sj.clone();
                extend_supports(&mut seg_j, &mut bookkeeper, w);

                let comp = if seg_j.logic.ends_with('*') { "".to_string() } else { "*".to_string() };
                let seg_i = Segment {
                    left_support: seg_j.right_support.iter().take(w).rev().map(|s| format!("{}{}", s.trim_end_matches('*'), comp)).collect(),
                    logic: format!("{}{}", seg_j.logic.trim_end_matches('*'), comp),
                    right_support: seg_j.left_support.iter().rev().take(w).map(|s| format!("{}{}", s.trim_end_matches('*'), comp)).collect(),
                    timer: format!("T{}", i),
                };
                seq[i] = Some(seg_i);
                seq[j] = Some(seg_j);
            }
            (Some(si), Some(sj)) => {
                let mut seg_i = si.clone();
                let mut seg_j = sj.clone();
                assert_eq!(seg_i.logic.trim_end_matches('*'), seg_j.logic.trim_end_matches('*'));
                assert!(seg_i.logic.ends_with('*') || seg_j.logic.ends_with('*'));
                assert!(!(seg_i.logic.ends_with('*') && seg_j.logic.ends_with('*')));
                if seg_i.left_support.len() < seg_j.left_support.len() {
                    assert!(seg_i.right_support.len() < seg_j.right_support.len());
                    // TODO: assert that existing support matches before we overwrite it.
                    extend_supports(&mut seg_j, &mut bookkeeper, w);
                    if seg_i.left_support.len() + seg_i.right_support.len() < w {
                        let comp = if seg_j.logic.ends_with('*') { "".to_string() } else { "*".to_string() };
                        seg_i = Segment {
                            left_support: seg_j.right_support.iter().take(w).rev().map(|s| format!("{}{}", s.trim_end_matches('*'), comp)).collect(),
                            logic: format!("{}{}", seg_j.logic.trim_end_matches('*'), comp),
                            right_support: seg_j.left_support.iter().rev().take(w).map(|s| format!("{}{}", s.trim_end_matches('*'), comp)).collect(),
                            timer: format!("T{}", i),
                        };
                    }
                } else if seg_j.left_support.len() < seg_i.left_support.len() {
                    assert!(seg_j.right_support.len() < seg_i.right_support.len());
                    // TODO: assert that existing support matches before we overwrite it.
                    extend_supports(&mut seg_i, &mut bookkeeper, w);
                    if seg_j.left_support.len() + seg_j.right_support.len() < w {
                        let comp = if seg_i.logic.ends_with('*') { "".to_string() } else { "*".to_string() };
                        seg_j = Segment {
                            left_support: seg_i.right_support.iter().take(w).rev().map(|s| format!("{}{}", s.trim_end_matches('*'), comp)).collect(),
                            logic: format!("{}{}", seg_i.logic.trim_end_matches('*'), comp),
                            right_support: seg_i.left_support.iter().rev().take(w).map(|s| format!("{}{}", s.trim_end_matches('*'), comp)).collect(),
                            timer: format!("T{}", j),
                        };
                    }
                }
                seq[i] = Some(seg_i);
                seq[j] = Some(seg_j);
            }
        }
    }

    seq
}



#[cfg(test)]
mod tests {
    use super::*;
    use crate::acfps::Acfp;
    use crate::checks::dlfolding::{build_pair_scores, nussinov, traceback_structures};
    use structure::{DotBracketVec, PairTable};

    #[test]
    fn test_assign_segments_balancing_and_expansion() {
        // Define the pair hierarchy
        let pair_hierarchy: FxHashMap<(usize, usize), usize> = vec![
            ((1, 2), 1),
            ((3, 4), 2),
            ((2, 3), 3),
        ].into_iter().collect();

        // One connected component
        let components = vec![vec![1, 2, 3, 4]];

        // Call the function
        let segments = assign_segments(&components, &pair_hierarchy);
        let sequence: Vec<String> = segments.iter().flatten().flat_map(|d| d.full_sequence()).collect();
        println!("{}", sequence.join(" "));

        // Print results for debugging
        for seg in segments.iter().flatten() {
            println!("{:?}", seg.full_sequence());
        }

        // Check symmetry and logic domain names
        let s1 = segments[1].as_ref().unwrap();
        let s2 = segments[2].as_ref().unwrap();
        let s3 = segments[3].as_ref().unwrap();
        let s4 = segments[4].as_ref().unwrap();

        // Ensure logic domains are complementary
        assert_eq!(s1.logic, "L1");
        assert_eq!(s2.logic, "L1*");
        assert_eq!(s3.logic, "L1");
        assert_eq!(s4.logic, "L1*");

        // Ensure left + right support lengths = 2w
        assert_eq!(s1.left_support.len() + s1.right_support.len(), 2);
        assert_eq!(s2.left_support.len() + s2.right_support.len(), 6);
        assert_eq!(s3.left_support.len() + s3.right_support.len(), 6);
        assert_eq!(s4.left_support.len() + s4.right_support.len(), 4);
    }

    #[test]
    fn test_assign_acfp_04() {
        let acfp = Acfp::try_from(". () (). ()() (().) ").unwrap();
        assert!(acfp.is_valid());

        let pairh = acfp.pair_hierarchy().unwrap();
        let ccomp = acfp.connected_components();
        println!("{:?}", pairh);
        println!("{:?}", ccomp);

        let pairh: FxHashMap<(usize, usize), usize> = vec![
            ((1, 2), 3),
            ((3, 4), 1),
            ((2, 3), 2),
            ((1, 5), 4),
        ].into_iter().collect();
        let ccomp = vec![vec![0], vec![1, 2, 3, 4, 5]];

        println!("{:?}", pairh);
        println!("{:?}", ccomp);

        let segments = assign_segments(&ccomp, &pairh);

        let mut transcript: Vec<String> = vec![];
        let mut tr_lengths: Vec<usize> = vec![];
        let mut ld_indices: Vec<usize> = vec![];
        for seg in segments.iter().flatten() {
            //println!("{:?}", seg.full_sequence());
            let ldi = transcript.len() + seg.left_support().len();
            transcript.extend(seg.full_sequence());
            tr_lengths.push(transcript.len()-1);
            ld_indices.push(ldi);
            println!("{:?} {} {}", transcript, transcript.len(), transcript[ldi]);
        }

        let sequence: Vec<String> = segments.iter().flatten().flat_map(|d| d.full_sequence()).collect();
        println!("{}", sequence.join(" "));
        println!("{}", transcript.join(" "));

        let p = build_pair_scores(&transcript, None);
        let dp = nussinov(&p);
        for l in tr_lengths {
            let all_structs = traceback_structures(0, l, &dp, &p);
            assert_eq!(all_structs.len(), 1);
            println!("{}", DotBracketVec::try_from(&PairTable(all_structs[0].clone())).unwrap());
        }

        assert!(false);
    }


}
