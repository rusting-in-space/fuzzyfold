use structure::PairTable;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Base { A, C, G, U }

impl TryFrom<char> for Base {
    type Error = ();
    fn try_from(c: char) -> Result<Self, ()> {
        Ok(match c.to_ascii_uppercase() {
            'A' => Base::A,
            'C' => Base::C,
            'G' => Base::G,
            'U' | 'T' => Base::U,
            _ => return Err(()),
        })
    }
}

pub fn basify(seq: &str) -> Vec<Base> {
    seq.chars().map(Base::try_from).collect::<Result<_, _>>().unwrap()
}

pub struct SecondaryStructure {
    pub sequence: Vec<Base>,
    pub pairings: PairTable,
}

impl SecondaryStructure {
    pub fn new(sequence: &str, pairings: &str) -> Result<Self, ()> {
        let seq = basify(sequence);
        let pt = PairTable::try_from(pairings).map_err(|_| ())?;
        Ok(SecondaryStructure { sequence: seq, pairings: pt })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_secondary_structure_new() {
        // Simple hairpin:  G C ( ( A A A ) ) G C
        let seq = "GCAAAGC";
        let dbn = "((...))"; // dot-bracket notation

        let ss = SecondaryStructure::new(seq, dbn).expect("should parse");
        assert_eq!(ss.sequence.len(), 7);
        assert_eq!(ss.sequence[0], Base::G);
        assert_eq!(ss.sequence[3], Base::A);
        assert_eq!(ss.pairings.len(), 7);
        assert_eq!(ss.pairings[0], Some(6));
        assert_eq!(ss.pairings[6], Some(0));
    }
}

 
