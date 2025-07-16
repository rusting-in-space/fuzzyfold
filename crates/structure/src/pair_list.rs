
use std::convert::TryFrom;
use crate::pair_table::PairTable;
use crate::error::StructureError;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Pair(pub usize, pub usize);

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairList {
    pub length: usize,
    pub pairs: Vec<Pair>,
}

impl TryFrom<&PairTable> for PairList {
    type Error = StructureError;

    fn try_from(pt: &PairTable) -> Result<Self, Self::Error> {
        match pt {
            PairTable::Single(pt) => {
                let mut pairs = Vec::new();
                for (i, &j_opt) in pt.iter().enumerate() {
                    if let Some(j) = j_opt {
                        if j == i {
                            return Err(StructureError::InvalidToken("self-pairing".to_string(), i));
                        }
                        if j >= pt.len() {
                            return Err(StructureError::InvalidToken(
                                    format!("invalid index {} at {}", j, i),
                                    i,
                            ));
                        }
                        match pt.get(j) {
                            Some(&Some(k)) if k == i => {
                                if i < j {
                                    pairs.push(Pair(i + 1, j + 1));
                                }
                            }
                            _ => return Err(StructureError::UnmatchedOpen(i)),
                        }
                    }
                }
                Ok(PairList {
                    length: pt.len(),
                    pairs,
                })
            }
            PairTable::Multi(_) => { 
                todo!();
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pair_list_from_pair_table() {
        let pt = PairTable::try_from("((..))").unwrap();
        let pl = PairList::try_from(&pt).unwrap();

        assert_eq!(pl.length, 6);
        assert_eq!(pl.pairs, vec![Pair(1, 6), Pair(2, 5)]); // 1-based
    }
}
