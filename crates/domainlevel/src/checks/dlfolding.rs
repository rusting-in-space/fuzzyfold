use ndarray::Array2;
use rustc_hash::{FxHashMap, FxHashSet};
use crate::utils::is_complement;

pub fn nussinov(p: &Array2<usize>) -> Array2<usize> {
    let (n, m) = p.dim();
    assert!(n == m);
    let mut dp = Array2::from_elem((n, n), 0);
    for l in 1..n {
        for i in 0..n - l {
            let j = i + l;
            let mut max_val = dp[(i + 1, j)].max(dp[(i, j - 1)]);
            if p[(i, j)] > 0 {
                max_val = max_val.max(dp[(i + 1, j - 1)] + p[(i, j)]);
            }
            for k in i + 1..j {
                max_val = max_val.max(dp[(i, k)] + dp[(k + 1, j)]);
            }
            dp[(i, j)] = max_val;
        }
    }
    dp
}

/// Returns a pairwise score matrix for a list of domain strings.
/// If a `lengths` dictionary is provided, it overrides the default score (1)
/// with the entry-specific length for each domain.
pub fn build_pair_scores(
    domains: &[String],
    lengths: Option<&FxHashMap<String, usize>>,
) -> Array2<usize> {
    let n = domains.len();
    let mut p = Array2::from_elem((n, n), 0);

    for ((i, j), value) in p.indexed_iter_mut() {
        assert_eq!(*value, 0); // sanity check

        if is_complement(&domains[i], &domains[j]) {
            *value = match lengths {
                Some(map) => {
                    let li = map.get(&domains[i]);
                    let lj = map.get(&domains[j]);

                    // Use min of lengths if both found, fallback to one, then default to 1
                    match (li, lj) {
                        (Some(&li), Some(&lj)) => li.min(lj),
                        (Some(&l), None) | (None, Some(&l)) => l,
                        (None, None) => 1,
                    }
                }
                None => 1,
            };
        }
    }

    p
}

#[allow(dead_code)] // for MFE prediction, no subopts.
fn traceback(
    i: usize,
    j: usize,
    dp: &Array2<usize>,
    p: &Array2<usize>,
    pairs: &mut Vec<(usize, usize)>,
) {
    if i >= j {
        return;
    }
    let dp_ij = dp[(i,j)];

    if dp_ij == dp[(i + 1, j)] {
        traceback(i + 1, j, dp, p, pairs);
    } else if dp_ij == dp[(i, j - 1)] {
        traceback(i, j - 1, dp, p, pairs);
    } else if p[(i, j)] > 0 && dp_ij == dp[(i + 1, j - 1)] + p[(i, j)] {
        pairs.push((i, j));
        traceback(i + 1, j - 1, dp, p, pairs);
    } else {
        for k in i + 1..j {
            if dp_ij == dp[(i, k)] + dp[(k + 1, j)] {
                traceback(i, k, dp, p, pairs);
                traceback(k + 1, j, dp, p, pairs);
                break;
            }
        }
    }
}

pub fn traceback_structures(
    i: usize,
    j: usize,
    dp: &Array2<usize>,
    p: &Array2<usize>,
) -> Vec<Vec<Option<usize>>> {
 
    let as_pairs = traceback_all(i, j, dp, p);

    as_pairs.into_iter()
        .map(|ps| {
            let mut table = vec![None; j + 1];
            for (p, q) in ps {
                table[p] = Some(q);
                table[q] = Some(p);
            }
            table
        })
        .collect()
}

pub fn traceback_all(
    i: usize,
    j: usize,
    dp: &Array2<usize>,
    p: &Array2<usize>,
) -> Vec<Vec<(usize, usize)>> {
    if i >= j {
        return vec![vec![]];
    }

    let mut results = vec![];
    let dp_ij = dp[(i,j)];

    // Case 1: i unpaired
    if dp_ij == dp[(i + 1, j)] {
        for sub in traceback_all(i + 1, j, dp, p) {
            results.push(sub);
        }
    }

    // Case 2: j unpaired
    if dp_ij == dp[(i, j - 1)] {
        for sub in traceback_all(i, j - 1, dp, p) {
            results.push(sub);
        }
    }

    // Case 3: i-j paired
    if p[(i, j)] > 0 && dp_ij == dp[(i + 1, j - 1)] + p[(i, j)] {
        for mut sub in traceback_all(i + 1, j - 1, dp, p) {
            sub.push((i, j));
            results.push(sub);
        }
    }

    // Case 4: bifurcation
    for k in i + 1..j {
        if dp_ij == dp[(i, k)] + dp[(k + 1, j)] {
            for left in traceback_all(i, k, dp, p) {
                for right in traceback_all(k + 1, j, dp, p) {
                    let mut combined = left.clone();
                    combined.extend(right);
                    results.push(combined);
                }
            }
        }
    }

    let mut seen: FxHashSet<Vec<(usize, usize)>> = FxHashSet::default();
    results
        .into_iter()
        .filter_map(|mut v| {
            v.sort();
            if seen.insert(v.clone()) {
                Some(v)
            } else {
                None
            }
        })
    .collect()
}


#[cfg(test)]
mod tests {
    use super::*;

    fn seq(domains: &[&str]) -> Vec<String> {
        domains.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn test_pair_score_simple() {
        let domains = seq(&["a", "a*", "b", "b*", "c"]);
        let p = build_pair_scores(&domains, None);
        assert_eq!(p[(0, 1)], 1);
        assert_eq!(p[(1, 0)], 1);
        assert_eq!(p[(2, 3)], 1);
        assert_eq!(p[(3, 2)], 1);
        assert_eq!(p[(0, 2)], 0);
    }

    #[test]
    fn test_nussinov_basic_structure() {
        let sequence = seq(&["a", "a*", "b", "b*"]);
        let sm: FxHashMap<String, usize> = [("a".to_string(), 2)].into_iter().collect();
        let p = build_pair_scores(&sequence, Some(&sm));
        let dp = nussinov(&p);
        assert_eq!(dp[(0, 3)], 3); // a-a* and b-b*

        let mut pairs: Vec<(usize, usize)> = Vec::new();
        traceback(0, sequence.len() - 1, &dp, &p, &mut pairs);
        assert_eq!(pairs, vec![(0, 1), (2, 3)]); // a-a* and b-b*
    }

    #[test]
    fn test_traceback_all_variants() {
        let sequence = seq(&["a", "x", "a*"]);
        let mut p = build_pair_scores(&sequence, None);
        p[(0, 2)] = 1; // ensure complement
        p[(2, 0)] = 1;
        let dp = nussinov(&p);
        let structs = traceback_all(0, sequence.len() - 1, &dp, &p);
        assert_eq!(dp[(0, 2)], 1);
        assert_eq!(structs.len(), 1);
        assert!(structs[0].contains(&(0, 2)));
    }

    #[test]
    fn test_traceback_all_bifurcation() {
        let sequence = seq(&["a", "a*", "a", "a*"]);
        let p = build_pair_scores(&sequence, None);
        let dp = nussinov(&p);
        let all_structs = traceback_all(0, sequence.len() - 1, &dp, &p);
        println!("{:?}", all_structs);
        assert_eq!(dp[(0, 3)], 2);

        assert_eq!(all_structs, [[(0, 3), (1, 2)], [(0, 1), (2, 3)]]);
        assert!(all_structs.iter().any(|s| s.contains(&(0, 1)) && s.contains(&(2, 3))));
        assert!(all_structs.iter().any(|s| s.contains(&(0, 3)) && s.contains(&(1, 2))));
        let all_structs = traceback_structures(0, sequence.len() - 1, &dp, &p);
        assert_eq!(all_structs, [[Some(3), Some(2), Some(1), Some(0)], 
                                 [Some(1), Some(0), Some(3), Some(2)]]);
    }

    #[test]
    fn test_traceback_all_multioutput() {
        let sequence = seq(&["a", "a*", "a", "a*", "a", "a*", "a", "a*"]);
        let p = build_pair_scores(&sequence, None);
        let dp = nussinov(&p);
        let all_structs = traceback_all(0, sequence.len() - 1, &dp, &p);
        println!("{:?}", all_structs);
        assert_eq!(all_structs.len(), 14);
        let all_structs = traceback_structures(0, sequence.len() - 1, &dp, &p);
        for s in &all_structs {
            println!("{:?}", s);
        }
        assert_eq!(all_structs.len(), 14);
    }

}

