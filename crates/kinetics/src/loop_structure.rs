use nohash_hasher::IntMap;
use nohash_hasher::IntSet;
use rustc_hash::FxHashMap;

use crate::NearestNeighborLoop;
use crate::LoopDecomposition;
use crate::EnergyModel;
use crate::utils::Base;

pub struct LoopCache {
    add_pairs: FxHashMap<NearestNeighborLoop, Vec<(NearestNeighborLoop, NearestNeighborLoop, i32)>>,
}

#[derive(Debug)]
pub struct LoopStructure {
    loop_list: IntMap<usize, (NearestNeighborLoop, i32)>,
    l_indices: IntSet<usize>,
    
    pair_list: Vec<(usize, usize)>, // 0-indexed pairs
 
    // sequence index to loop_list_index 
    // for unpaired: all i of same loop point to same loop
    // for pairs, i: outer loop, j: inner loop
    loop_lookup: IntMap<usize, usize>, 

    // Neighbor Caching --> easy to update locally, fast sampling access
    loop_neighbors: IntMap<usize, Vec<(usize, usize, i32)>>, // loop_list index to list of (i, j, deltaE)
    pair_neighbors: IntMap<usize, i32>, // pair_list index to deltaE
}

pub fn unpaired_pairs_in_loop<E: EnergyModel>(
    nn_loop: &NearestNeighborLoop,
    sequence: &[Base],
    model: &E,
) -> Vec<(usize, usize)> {

    let unpaired: Vec<usize> = match nn_loop {
        NearestNeighborLoop::Hairpin { closing: (i, j) } => {
            (*i + 1..*j).collect()
        }
        NearestNeighborLoop::Interior { closing: (i, j), inner: (k, l) } => {
            let mut v = Vec::new();
            v.extend(*i + 1..*k);
            v.extend(*l + 1..*j);
            v
        }
        NearestNeighborLoop::Multibranch { closing: (i, j), branches } => {
            let mut v = Vec::new();
            let mut start = *i;
            for &(p, q) in branches {
                v.extend(start + 1..p);
                start = q;
            }
            v.extend(start+1..=*j);
            v
        }
        NearestNeighborLoop::Exterior { branches } => {
            let mut v = Vec::new();
            let mut start = 0;
            v.push(start);
            for &(p, q) in branches {
                v.extend(start+1..p);
                start = q;
            }
            v.extend(start+1..sequence.len());
            v
        }
    };

    // Pair all (i, j) with j > i + 3 and check if pairable
    let mut pairs = Vec::new();
    for (idx_i, &i) in unpaired.iter().enumerate() {
        for &j in &unpaired[idx_i + 1..] {
            if j <= i + model.min_hp_size() {
                continue;
            }
            if model.can_pair(sequence[i], sequence[j]) {
                pairs.push((i, j));
            }
        }
    }

    pairs
}



impl LoopStructure {
    /// Produce a new updated SecStruct.
    pub fn apply_move(&mut self, i: usize, j: usize) {
        todo!("modify in place, depends on loop_neighbors and pair_neighbors.");
    }

    fn get_loop_neighbors<E: EnergyModel>(&mut self, cache: Option<LoopCache>, sequence: &[Base], model: &E) {
        for (idx, (nn_loop, en)) in self.loop_list.iter() {
            let pairs = unpaired_pairs_in_loop(nn_loop, sequence, model);
            println!("{:?}", pairs);
        }
    }
    

}


impl<T: LoopDecomposition, M: EnergyModel> TryFrom<(&[Base], &T, &M)> for LoopStructure {
    type Error = String;

    fn try_from((sequence, pairings, model): (&[Base], &T, &M)) -> Result<Self, Self::Error> {
        let mut loop_list = IntMap::default();
        let mut pair_list = Vec::new();
        let mut loop_lookup = IntMap::default();
        let mut idx = 0;

        pairings.for_each_loop(|l| {
            let energy = model.energy_of_loop(sequence, l);
            loop_list.insert(idx, (l.clone(), energy));

            match l {
                NearestNeighborLoop::Hairpin { closing: (i, j) } => {
                    pair_list.push((*i, *j));
                    for k in *i+1..=*j {
                        loop_lookup.insert(k, idx);
                    }
                }
                NearestNeighborLoop::Interior { closing: (i, j),  inner: (p, q) } => { 
                    pair_list.push((*i, *j));
                    for k in *i+1..=*p {
                        loop_lookup.insert(k, idx);
                    }
                    for k in *q+1..=*j {
                        loop_lookup.insert(k, idx);
                    }
                }
                NearestNeighborLoop::Multibranch { closing: (i, j), branches } => {
                    pair_list.push((*i, *j));
                    let mut start = *i;
                    for &(p, q) in branches {
                        for k in start+1..=p {
                            loop_lookup.insert(k, idx);
                        }
                        start = q;
                    }
                    for k in start+1..=*j {
                        loop_lookup.insert(k, idx);
                    }
                }
                NearestNeighborLoop::Exterior { branches } => {
                    let mut start = 0;
                    loop_lookup.insert(start, idx);
                    for &(p, q) in branches {
                        for k in start+1..=p {
                            loop_lookup.insert(k, idx);
                        }
                        start = q;
                    }
                    for k in start+1..sequence.len() {
                        loop_lookup.insert(k, idx);
                    }
                }
            }
            idx += 1;

        });

        Ok(LoopStructure {
            loop_list,
            l_indices: IntSet::default(),
            pair_list,
            loop_lookup,
            loop_neighbors: IntMap::default(),
            pair_neighbors: IntMap::default(),
        })
    }

}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::ViennaRNA;
    use structure::PairTable;

    #[test]
    fn test_loop_structure() {
        let sequence = "GCCCCGGUCA";
        let struct_0 = "...........";
        let struct_1 = ".(....)....";

        let model = ViennaRNA::new();
        let sequence: Vec<Base> = sequence.chars().map(Base::try_from).collect::<Result<_, _>>().unwrap();
        let struct_0 = PairTable::try_from(struct_0).expect("valid");
        let ls = LoopStructure::try_from((&sequence[..], &struct_0, &model)).expect("bla");
        println!("{:?}", ls);


        let struct_1 = PairTable::try_from(struct_1).expect("valid");
        let mut ls = LoopStructure::try_from((&sequence[..], &struct_0, &model)).expect("bla");
        println!("{:?}", ls);
        ls.get_loop_neighbors(None, &sequence, &model);

        assert!(false);
    }
}

