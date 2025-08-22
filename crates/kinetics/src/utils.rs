use structure::PairTable;
use crate::energy_tables::Base;
use crate::energy_tables::basify;

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

 
