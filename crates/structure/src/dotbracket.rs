use std::fmt;
use std::ops::Deref;
use std::convert::TryFrom;

use crate::pair_table::PairTable;
use crate::error::StructureError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DotBracket {
    Unpaired, // '.'
    Open,     // '('
    Close,    // ')'
}

impl From<char> for DotBracket {
    fn from(c: char) -> Self {
        match c {
            '.' => DotBracket::Unpaired,
            '(' => DotBracket::Open,
            ')' => DotBracket::Close,
            _ => panic!("Invalid dot-bracket character: {}", c),
        }
    }
}

impl From<DotBracket> for char {
    fn from(db: DotBracket) -> Self {
        match db {
            DotBracket::Open => '(',
            DotBracket::Close => ')',
            DotBracket::Unpaired => '.',
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DotBracketVec(pub Vec<DotBracket>);

impl Deref for DotBracketVec {
    type Target = [DotBracket];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<&str> for DotBracketVec {
    fn from(s: &str) -> Self {
        let vec = s.chars().map(DotBracket::from).collect();
        DotBracketVec(vec)
    }
}

impl TryFrom<&PairTable> for DotBracketVec {
    type Error = StructureError;

    fn try_from(pt: &PairTable) -> Result<Self, Self::Error> {
        let mut result = Vec::with_capacity(pt.len());

        for (i, &j_opt) in pt.iter().enumerate() {
            match j_opt {
                None => result.push(DotBracket::Unpaired),
                Some(j) => {
                    if j > i {
                        match pt.get(j) {
                            Some(Some(k)) if *k == i => {
                                result.push(DotBracket::Open);
                            }
                            _ => return Err(StructureError::UnmatchedOpen(i)),
                        }
                    } else if j < i {
                        match pt.get(j) {
                            Some(Some(k)) if *k == i => {
                                result.push(DotBracket::Close);
                            }
                            _ => return Err(StructureError::UnmatchedClose(i)),
                        }
                    } else { // j == i
                        return Err(StructureError::InvalidToken("self-pairing".into(), i));
                    }
                }
            }
        }

        Ok(DotBracketVec(result))
    }
}

impl fmt::Display for DotBracketVec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for db in &self.0 {
            write!(f, "{}", char::from(*db))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dot_bracket_from_char() {
        assert_eq!(DotBracket::from('.'), DotBracket::Unpaired);
        assert_eq!(DotBracket::from('('), DotBracket::Open);
        assert_eq!(DotBracket::from(')'), DotBracket::Close);
    }

    #[test]
    fn test_char_from_dot_bracket() {
        assert_eq!(char::from(DotBracket::Unpaired), '.');
        assert_eq!(char::from(DotBracket::Open), '(');
        assert_eq!(char::from(DotBracket::Close), ')');
    }

    #[test]
    #[should_panic(expected = "Invalid dot-bracket character")]
    fn test_dot_bracket_from_invalid_char() {
        let _ = DotBracket::from('x');
    }

    #[test]
    fn test_dot_bracket_vec_from_str() {
        let dbv = DotBracketVec::from("(.).");
        assert_eq!(format!("{}", dbv), "(.).");
        assert_eq!(dbv.len(), 4);
        assert_eq!(dbv[0], DotBracket::Open);
        assert_eq!(dbv[1], DotBracket::Unpaired);
        assert_eq!(dbv[2], DotBracket::Close);
        assert_eq!(dbv[3], DotBracket::Unpaired);
    }

    #[test]
    fn test_dot_bracket_vec_from_pair_table() {
        let pt = PairTable::try_from("((..))").unwrap();
        let dbv = DotBracketVec::try_from(&pt).unwrap();
        assert_eq!(format!("{}", dbv), "((..))");
    }

}

