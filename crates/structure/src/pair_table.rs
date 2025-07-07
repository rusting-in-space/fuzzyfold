use std::ops::{Deref, DerefMut};
use std::convert::TryFrom;
use crate::error::StructureError;
use crate::dotbracket::{DotBracket, DotBracketVec};
use crate::pair_list::{Pair, PairList};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairTable(pub Vec<Option<usize>>);

impl Deref for PairTable {
    type Target = [Option<usize>];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for PairTable {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}


impl TryFrom<&str> for PairTable {
    type Error = StructureError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut stack = Vec::new();
        let mut table = vec![None; s.len()];

        for (i, c) in s.chars().enumerate() {
            match c {
                '(' => stack.push(i),
                ')' => {
                    let j = stack.pop().ok_or(StructureError::UnmatchedClose(i))?;
                    table[i] = Some(j);
                    table[j] = Some(i);
                }
                '.' => (),
                other => return Err(StructureError::InvalidToken(other.to_string(), i)),
            }
        }

        if let Some(i) = stack.pop() {
            return Err(StructureError::UnmatchedOpen(i));
        }

        Ok(PairTable(table))
    }
}

impl TryFrom<&DotBracketVec> for PairTable {
    type Error = StructureError;

    fn try_from(db: &DotBracketVec) -> Result<Self, Self::Error> {
        let mut stack: Vec<usize> = Vec::new();
        let mut table = vec![None; db.len()];

        for (i, dot) in db.iter().enumerate() {
            match dot {
                DotBracket::Open => stack.push(i),
                DotBracket::Close => {
                    let j = stack.pop().ok_or(StructureError::UnmatchedClose(i))?;
                    if i == j {
                        return Err(StructureError::InvalidToken("self_pairing".to_string(), i));
                    }
                    table[i] = Some(j);
                    table[j] = Some(i);
                }
                DotBracket::Unpaired => {} // leave as None
            }
        }

        if let Some(i) = stack.pop() {
            return Err(StructureError::UnmatchedOpen(i));
        }

        Ok(PairTable(table))
    }
}


impl TryFrom<&PairList> for PairTable {
    type Error = StructureError;

    fn try_from(pl: &PairList) -> Result<Self, Self::Error> {
        let mut table = vec![None; pl.length];

        for &Pair(i1, j1) in &pl.pairs {
            // Convert 1-based to 0-based
            if i1 == j1 {
                return Err(StructureError::InvalidToken("self-pairing".to_string(), i1));
            }
            if i1 == 0 || j1 == 0 || i1 > pl.length || j1 > pl.length {
                return Err(StructureError::InvalidToken(
                    format!("Invalid 1-based pair ({}, {})", i1, j1),
                    i1.max(j1).saturating_sub(1),
                ));
            }

            let (i, j) = (i1 - 1, j1 - 1);

            if let Some(prev) = table[i] {
                if prev != j {
                    return Err(StructureError::InvalidToken(
                        format!("Conflict: {} already paired with {}", i + 1, prev + 1),
                        i,
                    ));
                }
            }
            if let Some(prev) = table[j] {
                if prev != i {
                    return Err(StructureError::InvalidToken(
                        format!("Conflict: {} already paired with {}", j + 1, prev + 1),
                        j,
                    ));
                }
            }

            table[i] = Some(j);
            table[j] = Some(i);
        }

        Ok(PairTable(table))
    }
}

impl PairTable {
    /// Check if the substructure from `i..j` is well-formed:
    /// - All pairings are internal to the interval
    pub fn is_well_formed(&self, i: usize, j: usize) -> bool {
        assert!(j <= self.len(), "Invalid interval: j must be <= length");

        for k in i..j {
            if let Some(l) = self[k] {
                if l < i || l >= j {
                    return false; // points outside
                }
            }
        }
        true
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_pair_table() {
        let pt = PairTable::try_from("((..))").unwrap();
        assert_eq!(pt.len(), 6);
        assert_eq!(pt[0], Some(5));
        assert_eq!(pt[1], Some(4));
        assert_eq!(pt[2], None);
        assert_eq!(pt[3], None);
        assert_eq!(pt[4], Some(1));
        assert_eq!(pt[5], Some(0));
    }

    #[test]
    fn test_unmatched_open() {
        let err = PairTable::try_from("(()").unwrap_err();
        assert_eq!(format!("{}", err), "Unmatched '(' at position 0");
    }

    #[test]
    fn test_unmatched_close() {
        let err = PairTable::try_from("())").unwrap_err();
        assert_eq!(format!("{}", err), "Unmatched ')' at position 2");
    }

    #[test]
    fn test_invalid_token() {
        let err = PairTable::try_from("(x)").unwrap_err();
        assert_eq!(format!("{}", err), "Invalid token 'x' at position 1");
    }

    #[test]
    fn test_well_formed_empty_interval() {
        let pt= PairTable::try_from("...").unwrap();
        assert!(pt.is_well_formed(0, 0)); 
        assert!(pt.is_well_formed(0, 1)); 
        assert!(pt.is_well_formed(0, 2)); 
        assert!(pt.is_well_formed(0, 3)); 
        assert!(pt.is_well_formed(1, 3)); 
        assert!(pt.is_well_formed(2, 3)); 
        assert!(pt.is_well_formed(3, 3)); 
    }

    #[test]
    fn test_well_formed_pairings_within_interval() {
        let pt = PairTable::try_from(".(.).").unwrap();
        assert!(pt.is_well_formed(0, 5)); // Full interval -- 0-based
        assert!(pt.is_well_formed(0, 4)); 
        assert!(pt.is_well_formed(1, 5));
        assert!(pt.is_well_formed(1, 4));
        assert!(pt.is_well_formed(1, 4));
        assert!(pt.is_well_formed(2, 3));
        assert!(!pt.is_well_formed(0, 3)); 
        assert!(!pt.is_well_formed(1, 3)); 
        assert!(!pt.is_well_formed(2, 4)); 
    }

    #[test]
    #[should_panic(expected = "Invalid interval: j must be <= length")]
    fn test_well_formed_out_of_bounds_assert() {
        let pt = PairTable::try_from("..").unwrap();
        pt.is_well_formed(0, 3); // j = pt.len(), should panic
    }

    #[test]
    fn test_pair_table_from_pair_list_valid() {
        let pl = PairList {
            length: 6,
            pairs: vec![Pair(1, 6), Pair(2, 5)],
        };
        let pt = PairTable::try_from(&pl).unwrap();
        assert_eq!(pt.0, vec![Some(5), Some(4), None, None, Some(1), Some(0)]);
    }

    #[test]
    fn test_pair_table_from_pair_list_out_of_bounds() {
        let pl = PairList {
            length: 4,
            pairs: vec![Pair(1, 5)],
        };
        let err = PairTable::try_from(&pl).unwrap_err();
        assert!(matches!(err, StructureError::InvalidToken(_, _)));
    }

    #[test]
    fn test_pair_table_from_pair_list_self_pairing() {
        let pl = PairList {
            length: 4,
            pairs: vec![Pair(2, 2)],
        };
        let err = PairTable::try_from(&pl).unwrap_err();
        assert!(matches!(err, StructureError::InvalidToken(..)));
    }

    #[test]
    fn test_pair_table_from_pair_list_conflicting() {
        let pl = PairList {
            length: 6,
            pairs: vec![Pair(1, 6), Pair(1, 5)],
        };
        let err = PairTable::try_from(&pl).unwrap_err();
        assert!(matches!(err, StructureError::InvalidToken(_, _)));
    }
}


