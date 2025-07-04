use std::fmt;

#[derive(Debug)]
pub enum StructureError {
    UnmatchedOpen(usize),          // '(' at this position was never closed
    UnmatchedClose(usize),         // ')' at this position has no matching '('
    InvalidToken(String, usize),   // invalid char and position
}

impl fmt::Display for StructureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            StructureError::UnmatchedOpen(i) => {
                write!(f, "Unmatched '(' at position {}", i)
            }
            StructureError::UnmatchedClose(i) => {
                write!(f, "Unmatched ')' at position {}", i)
            }
            StructureError::InvalidToken(tok, i) => {
                write!(f, "Invalid token '{}' at position {}", tok, i)
            }
        }
    }
}

impl std::error::Error for StructureError {}


