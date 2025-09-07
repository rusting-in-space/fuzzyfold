use nohash_hasher::IntMap;
use nohash_hasher::IntSet;
use rustc_hash::FxHashMap;

use energy::NearestNeighborLoop;
use energy::LoopDecomposition;
use energy::EnergyModel;
use energy::Base;

pub struct LoopCache {
    eval_loop: FxHashMap<NearestNeighborLoop, i32>,
    add_pairs: FxHashMap<NearestNeighborLoop, Vec<(NearestNeighborLoop, NearestNeighborLoop, i32)>>,
}

#[derive(Debug)]
pub struct LoopStructure<'a, M: EnergyModel> {
    sequence: &'a [Base],
    model: &'a M,

    loop_list: IntMap<usize, (NearestNeighborLoop, i32)>,
    l_indices: IntSet<usize>,
    
    pair_list: IntMap<usize, (usize, usize)>,
    p_indices: IntSet<usize>,
 
    // sequence index to loop_list_index 
    // for unpaired: all i of same loop point to same loop
    // for pairs, i: outer loop, j: inner loop
    loop_lookup: IntMap<usize, usize>, 

    // Neighbor Caching --> easy to update locally, fast sampling access
    loop_neighbors: IntMap<usize, Vec<(usize, usize, i32)>>, // loop_list index to list of (i, j, deltaE)
    pair_neighbors: IntMap<usize, i32>, // pair_list index to deltaE
}

fn allocate_index(indices: &mut IntSet<usize>, fallback: usize) -> usize {
    indices.iter().next().cloned().map(|x| {
        indices.remove(&x);
        x
    }).unwrap_or(fallback)
}


impl<'a, M: EnergyModel> LoopStructure<'a, M> {
    pub fn get_add_neighbors(&self) -> Vec<(usize, usize, i32)> {
        self.loop_neighbors
            .values()
            .flat_map(|v| v.iter().cloned())
            .collect()
    }

    pub fn get_del_neighbors(&self) -> Vec<(usize, (usize, usize), i32)> {
        self.pair_neighbors
            .iter()
            .map(|(pli, &delta_e)| (*pli, self.pair_list[pli], delta_e))
            .collect()
    }

    pub fn energy(&self) -> i32 {
        self.loop_list
            .values()
            .map(|(_, e)| *e)
            .sum()
    }

     /// Produce a new updated SecStruct.
    pub fn apply_del_move(&mut self, pl_idx: usize) {
        let &(i, j) = self.pair_list.get(&pl_idx).expect("Missing pair_list entry.");

        let &lli = self.loop_lookup.get(&i).expect("Missing loop_lookup entry for i.");
        let &llj = self.loop_lookup.get(&j).expect("Missing loop_lookup entry for j.");

        let (outer, o_en) = self.loop_list.get(&lli).expect("Missing loop_list entry lli.");
        let (inner, i_en) = self.loop_list.get(&llj).expect("Missing loop_list entry llj.");

        // Recompute or get from global cache.
        let combo = &outer.join_loop(&inner);
        let delta = self.pair_neighbors.get(&pl_idx).expect("Missing pair_neighbors entry.");
        let combo_energy = (o_en + i_en) - delta;

        // Cleanup and re-allocate.
        self.pair_list.remove(&pl_idx);
        self.pair_neighbors.remove(&pl_idx);
        self.p_indices.insert(pl_idx);

        self.loop_list.remove(&lli);
        self.loop_list.remove(&llj);
        self.loop_neighbors.remove(&lli);
        self.loop_neighbors.remove(&llj);
        self.l_indices.insert(llj);

        self.loop_list.insert(lli, (combo.clone(), combo_energy));
        let neighbors = get_loop_neighbors(&combo, &combo_energy, self.sequence, self.model);
        self.loop_neighbors.insert(lli, neighbors);

        for k in &combo.inclusive_unpaired_indices(self.sequence.len()) {
            assert!(self.loop_lookup[k] == lli || self.loop_lookup[k] == llj);
            self.loop_lookup.insert(*k, lli);
        }
    }

    /// Produce a new updated SecStruct.
    pub fn apply_add_move(&mut self, i: usize, j: usize) {
        let model = self.model;
        let &lli = self.loop_lookup.get(&i).expect("Missing loop_lookup entry for i.");

        assert_eq!(
            lli, 
            *self.loop_lookup.get(&j).expect("Missing loop_lookup entry for j.")
        );

        let (combo, c_energy) = self.loop_list.get(&lli).expect("Missing loop_list entry lli.");
        let (outer, inner) = combo.split_loop(i, j);
        //NOTE; could look that up.
        let outer_energy = model.energy_of_loop(self.sequence, &outer);
        let inner_energy = model.energy_of_loop(self.sequence, &inner);
        let delta = (outer_energy + inner_energy) - c_energy;
  
        self.loop_list.remove(&lli);
        self.loop_neighbors.remove(&lli);
        self.loop_list.insert(lli, (outer.clone(), outer_energy));
        self.loop_neighbors.insert(lli, 
            get_loop_neighbors(&outer, &outer_energy, self.sequence, model));

        let llj = allocate_index(&mut self.l_indices, self.loop_list.len());

        self.loop_list.insert(llj, (inner.clone(), inner_energy));
        self.loop_neighbors.insert(llj, 
            get_loop_neighbors(&inner, &inner_energy, self.sequence, model));

        let pli = allocate_index(&mut self.p_indices, self.pair_list.len());
        self.pair_list.insert(pli, (i, j));

        for k in &outer.inclusive_unpaired_indices(self.sequence.len()) {
            self.loop_lookup.insert(*k, lli);
        }

        for k in &inner.inclusive_unpaired_indices(self.sequence.len()) {
            self.loop_lookup.insert(*k, llj);
        }

        self.pair_neighbors.insert(pli, delta);

    }

}

impl<'a, T: LoopDecomposition, M: EnergyModel> TryFrom<(&'a [Base], &T, &'a M)> for LoopStructure<'a, M> {
    type Error = String;

    fn try_from((sequence, pairings, model): (&'a [Base], &T, &'a M)) -> Result<Self, Self::Error> {
        let mut loop_list = IntMap::default();
        let mut pair_list = IntMap::default();
        let mut loop_lookup = IntMap::default();
        let mut ll_idx = 0;
        let mut pl_idx = 0;

        pairings.for_each_loop(|l| {
            let energy = model.energy_of_loop(sequence, l);
            loop_list.insert(ll_idx, (l.clone(), energy));
            if let Some((i, j)) = l.closing() {
                pair_list.insert(pl_idx, (i, j)); 
                pl_idx += 1;
            }
            for k in &l.inclusive_unpaired_indices(sequence.len()) {
                loop_lookup.insert(*k, ll_idx);
            }
            ll_idx += 1;
        });

        let mut loop_neighbors = IntMap::default();
        for (nn_idx, (nn_loop, nn_energy)) in loop_list.iter() {
            let neighbors = get_loop_neighbors(&nn_loop, nn_energy, sequence, model);
            loop_neighbors.insert(*nn_idx, neighbors);
        }

        let mut pair_neighbors = IntMap::default();
        for (pl_idx, (i, j)) in pair_list.iter() {
            let (outer, outer_energy) = loop_list.get(loop_lookup.get(i).unwrap()).unwrap();
            let (inner, inner_energy) = loop_list.get(loop_lookup.get(j).unwrap()).unwrap();
            let combo = &outer.join_loop(&inner);
            let combo_energy = model.energy_of_loop(sequence, &combo);
            let delta = combo_energy - (outer_energy + inner_energy);
            pair_neighbors.insert(*pl_idx, delta);
        }

        Ok(LoopStructure {
            sequence,
            model,
            loop_list,
            l_indices: IntSet::default(),
            pair_list,
            p_indices: IntSet::default(),
            loop_lookup,
            loop_neighbors,
            pair_neighbors,
        })
    }

}

fn get_loop_neighbors<E: EnergyModel>(
    nn_loop: &NearestNeighborLoop, 
    nn_energy: &i32, 
    sequence: &[Base], 
    model: &E
) -> Vec<(usize, usize, i32)> {
    let mut result = vec![];
    for (i, j) in unpaired_pairs_in_loop(nn_loop, sequence, model) {
        let (outer, inner) = nn_loop.split_loop(i, j);
        let outer_energy = model.energy_of_loop(sequence, &outer);
        let inner_energy = model.energy_of_loop(sequence, &inner);
        let delta = (outer_energy + inner_energy) - nn_energy;
        result.push((i, j, delta));
    }
    result
}

pub fn unpaired_pairs_in_loop<E: EnergyModel>(
    nn_loop: &NearestNeighborLoop,
    sequence: &[Base],
    model: &E,
) -> Vec<(usize, usize)> {
    let unpaired = nn_loop.unpaired_indices(sequence.len());
    // Pair all (i, j) with j > i + 3 and check if pairable
    let mut pairs = Vec::new();
    for (idx_i, &i) in unpaired.iter().enumerate() {
        for &j in &unpaired[idx_i + 1..] {
            if j <= i + model.min_hairpin_size() {
                continue;
            }
            if model.can_pair(sequence[i], sequence[j]) {
                pairs.push((i, j));
            }
        }
    }
    pairs
}



#[cfg(test)]
mod tests {
    use super::*;
    use structure::PairTable;
    use energy::ViennaRNA;
    use energy::basify;

    #[test]
    fn test_loop_structure() {
        let sequence = "GCCCCGGUCA";
        let struct_0 = "...........";
        let struct_1 = ".(....)....";

        let model = ViennaRNA::default();
        let sequence = basify(sequence);
        let struct_0 = PairTable::try_from(struct_0).expect("valid");
        let mut ls = LoopStructure::try_from((&sequence[..], &struct_0, &model)).expect("bla");

        println!("{:?} {}", ls, ls.energy());
        println!("Add neighbors");
        for (i, j, de) in ls.get_add_neighbors() {
            println!("{} {} {}", i, j, de);
        }

        println!("Del neighbors");
        for (pi, (i, j), de) in ls.get_del_neighbors() {
            println!("{} ({}, {}), {}", pi, i, j, de);
        }

        println!("Adddig 1, 5");
        ls.apply_add_move(1, 5);
        println!("{:?} {}", ls, ls.energy());
        println!("Add neighbors");
        for (i, j, de) in ls.get_add_neighbors() {
            println!("{} {} {}", i, j, de);
        }

        println!("Del neighbors");
        for (pi, (i, j), de) in ls.get_del_neighbors() {
            println!("{} ({}, {}), {}", pi, i, j, de);
        }

        println!("Removing 1, 5");
        ls.apply_del_move(0);
        println!("{:?} {}", ls, ls.energy());
        println!("Add neighbors");
        for (i, j, de) in ls.get_add_neighbors() {
            println!("{} {} {}", i, j, de);
        }

        println!("Del neighbors");
        for (pi, (i, j), de) in ls.get_del_neighbors() {
            println!("{} ({}, {}), {}", pi, i, j, de);
        }


        let struct_1 = PairTable::try_from(struct_1).expect("valid");
        let mut ls = LoopStructure::try_from((&sequence[..], &struct_1, &model)).expect("bla");
        println!("Add neighbors");
        for (i, j, de) in ls.get_add_neighbors() {
            println!("{} {} {}", i, j, de);
        }

        println!("Del neighbors");
        for (pi, (i, j), de) in ls.get_del_neighbors() {
            println!("{} ({}, {}), {}", pi, i, j, de);
        }


        //println!("{:?} {}", ls, ls.energy());

        //assert!(false);
    }
}

