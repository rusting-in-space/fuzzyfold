
use crate::NearestNeighborLoop;
use crate::LoopDecomposition;
use crate::energy_tables::Base;

use structure::PairTable;

pub trait EnergyModel {
    fn can_pair(&self, b1: Base, b2: Base) -> bool;

    fn min_hairpin_size(&self) -> usize;

    fn energy_of_structure<T: LoopDecomposition>(&self, 
        sequence: &[Base], 
        structure: &T
    ) -> i32;

    /// The free energy difference between forming the pair vs breaking it.
    //TODO: may not be particularly efficient, and depends on PairTable!
    fn energy_of_pair(&self, 
        sequence: &[Base], 
        structure: &PairTable,
        i: usize,
        j: usize,
    ) -> i32;

    fn energy_of_loop(&self, 
        sequence: &[Base], 
        nn_loop: &NearestNeighborLoop
    ) -> i32;
}

#[cfg(test)]
mod tests {
    use super::*;
    use Base::*;

    struct MockEnergyModel;

    impl EnergyModel for MockEnergyModel {
        fn can_pair(&self, b1: Base, b2: Base) -> bool {
            matches!((b1, b2), (A, U) | (U, A) | (C, G) | (G, C))
        }

        fn min_hairpin_size(&self) -> usize {
            3
        }

        fn energy_of_structure<T: LoopDecomposition>(
            &self,
            _sequence: &[Base],
            _structure: &T,
        ) -> i32 {
            -10 // dummy value
        }

        fn energy_of_pair(
            &self,
            _sequence: &[Base],
            _structure: &PairTable,
            _i: usize,
            _j: usize,
        ) -> i32 {
            -3 // dummy value
        }

        fn energy_of_loop(
            &self,
            _sequence: &[Base],
            _nn_loop: &NearestNeighborLoop,
        ) -> i32 {
            5 // dummy value
        }
    }

    #[test]
    fn test_can_pair() {
        let model = MockEnergyModel;

        assert!(model.can_pair(A, U));
        assert!(model.can_pair(C, G));
        assert!(!model.can_pair(A, G));
        assert!(!model.can_pair(C, C));
    }

    #[test]
    fn test_min_hairpin_size() {
        let model = MockEnergyModel;
        assert_eq!(model.min_hairpin_size(), 3);
    }

    #[test]
    fn test_energy_of_structure() {
        let model = MockEnergyModel;

        struct DummyStructure;
        impl LoopDecomposition for DummyStructure {
            fn for_each_loop<F: FnMut(&NearestNeighborLoop)>(&self, f: F) {

            }
            fn loop_enclosed_by(&self, closing: Option<(usize, usize)>) -> NearestNeighborLoop {
                NearestNeighborLoop::Hairpin { closing: closing.unwrap() }
            }
            fn get_enclosing_pair(&self, i: usize, j: usize) -> Option<(usize, usize)> {
                Some((i, j))
            }
        }

        let sequence = vec![A, U, G, C, U];
        let structure = DummyStructure;

        let energy = model.energy_of_structure(&sequence, &structure);
        assert_eq!(energy, -10);
    }

    #[test]
    fn test_energy_of_loop() {
        let model = MockEnergyModel;

        let sequence = vec![A, U, C, G];
        let nn_loop = NearestNeighborLoop::Hairpin { closing: (0, 3) };

        let energy = model.energy_of_loop(&sequence, &nn_loop);
        assert_eq!(energy, 5);
    }
}


