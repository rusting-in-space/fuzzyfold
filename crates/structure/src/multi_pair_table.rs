use std::convert::TryFrom;
use crate::StructureError;
use crate::DotBracket;
use crate::DotBracketVec;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultiPairTable(pub Vec<Vec<Option<(usize, usize)>>>);

impl MultiPairTable {
    pub fn len(&self) -> usize {
        self.0.iter().map(|s| s.len()).sum()
    }
    
    pub fn is_well_formed(&self, _loc1: (usize, usize), _loc2: (usize, usize)) -> bool {
        todo!("not implemented");
    }

    pub fn get_pair(&self, loc: (usize, usize)) -> &Option<(usize, usize)> {
        &self.0[loc.0][loc.1]
    }
}

impl TryFrom<&str> for MultiPairTable {
    type Error = StructureError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut strand_idx = 0;
        let mut domain_idx = 0;
        let mut stack: Vec<(usize, usize)> = Vec::new(); // (strand_idx, domain_idx)
        let mut pair_table: Vec<Vec<Option<(usize, usize)>>> = vec![vec![]];

        for (i, ch) in s.chars().enumerate() {
            match ch {
                '+' | '&' => {
                    assert!(domain_idx + strand_idx > 0);
                    if i < s.len()-1 {
                        // Don't append an empty vec in the "hack version"
                        pair_table.push(vec![]);
                    }
                    strand_idx += 1;
                    domain_idx = 0;
                }
                '(' => {
                    stack.push((strand_idx, domain_idx));
                    pair_table[strand_idx].push(None); // placeholder
                    domain_idx += 1;
                }
                ')' => {
                    let (si, di) = stack.pop().ok_or(StructureError::UnmatchedMultiClose((strand_idx, domain_idx)))?;
                    pair_table[si][di] = Some((strand_idx, domain_idx));
                    pair_table[strand_idx].push(Some((si, di)));
                    domain_idx += 1;
                }
                '.' => {
                    pair_table[strand_idx].push(None);
                    domain_idx += 1;
                }
                _ => {
                    return Err(StructureError::InvalidToken(ch.to_string(), i)); 
                }
            }
        }
        if let Some((si, di)) = stack.pop() {
            return Err(StructureError::UnmatchedMultiOpen((si, di)));
        }
        Ok(MultiPairTable(pair_table))
    }
}

impl TryFrom<&DotBracketVec> for MultiPairTable {
    type Error = StructureError;

    fn try_from(db: &DotBracketVec) -> Result<Self, Self::Error> {
        // Multi-stranded case
        let mut strand_idx = 0;
        let mut domain_idx = 0;
        let mut stack: Vec<(usize, usize)> = Vec::new(); // (strand_idx, domain_idx)
        let mut pair_table: Vec<Vec<Option<(usize, usize)>>> = vec![vec![]];

        for dot in db.iter() {
            match dot {
                DotBracket::Break => {
                    assert!(domain_idx + strand_idx > 0);
                    pair_table.push(vec![]);
                    strand_idx += 1;
                    domain_idx = 0;
                }
                DotBracket::Open => {
                    stack.push((strand_idx, domain_idx));
                    pair_table[strand_idx].push(None); // placeholder
                    domain_idx += 1;
                }
                DotBracket::Close => {
                    let (si, di) = stack.pop().ok_or(StructureError::UnmatchedMultiClose((strand_idx, domain_idx)))?;
                    pair_table[si][di] = Some((strand_idx, domain_idx));
                    pair_table[strand_idx].push(Some((si, di)));
                    domain_idx += 1;
                }
                DotBracket::Unpaired => {
                    pair_table[strand_idx].push(None);
                    domain_idx += 1;
                }
            }
        }

        if let Some((si, di)) = stack.pop() {
            return Err(StructureError::UnmatchedMultiOpen((si, di)));
        }

        Ok(MultiPairTable(pair_table))
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multi_pair_table() {
        let pt = MultiPairTable::try_from("((.+.))").unwrap();
        assert_eq!(pt.len(), 6);
        assert_eq!(*pt.get_pair((0, 0)), Some((1, 2)));
        assert_eq!(*pt.get_pair((0, 1)), Some((1, 1)));
        assert_eq!(*pt.get_pair((0, 2)), None);
        assert_eq!(*pt.get_pair((1, 0)), None);
        assert_eq!(*pt.get_pair((1, 1)), Some((0, 1)));
        assert_eq!(*pt.get_pair((1, 2)), Some((0, 0)));
    }

    #[test]
    fn test_multi_pair_table_hack() {
        let pt = MultiPairTable::try_from("((..))+").unwrap();
        assert_eq!(pt.len(), 6);
        assert_eq!(*pt.get_pair((0, 0)), Some((0, 5)));
        assert_eq!(*pt.get_pair((0, 1)), Some((0, 4)));
        assert_eq!(*pt.get_pair((0, 2)), None);
        assert_eq!(*pt.get_pair((0, 3)), None);
        assert_eq!(*pt.get_pair((0, 4)), Some((0, 1)));
        assert_eq!(*pt.get_pair((0, 5)), Some((0, 0)));
    }
}


