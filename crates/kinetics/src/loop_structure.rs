use std::fmt;
use nohash_hasher::IntMap;
use nohash_hasher::IntSet;
use rustc_hash::FxHashMap;

use structure::DotBracket;
use structure::DotBracketVec;
use energy::NearestNeighborLoop;
use energy::LoopDecomposition;
use energy::EnergyModel;
use energy::Base;

#[allow(dead_code)]
pub struct LoopCache {
    eval_loop: FxHashMap<NearestNeighborLoop, i32>,
    add_pairs: FxHashMap<NearestNeighborLoop, Vec<(NearestNeighborLoop, NearestNeighborLoop, i32)>>,
}

pub struct LoopStructure<'a, M: EnergyModel> {
    sequence: &'a [Base],
    model: &'a M,

    // They should be thought of "lists", but for easier updates,
    // they are actually implemented ad IntMaps with some index bookkeeping.
    // After all, we don't want to update the whole loop_lookup table after
    // modifying these "lists".
    loop_list: IntMap<usize, (NearestNeighborLoop, i32)>,
    l_indices: IntSet<usize>,
    
    pair_list: IntMap<usize, usize>,
 
    // sequence index to loop_list index 
    // for unpaired: all i of same loop point to same loop
    // for pairs, i: outer loop, j: inner loop
    loop_lookup: IntMap<usize, usize>, 

    // Neighbor Caching --> easy to update locally, fast sampling access
    loop_neighbors: IntMap<usize, Vec<(usize, usize, i32)>>, // loop_list index to list of (i, j, deltaE)
    pair_neighbors: IntMap<usize, i32>, // pair_list index to deltaE
}

impl<'a, M: EnergyModel> fmt::Debug for LoopStructure<'a, M> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LoopStructure")
            .field("num_loops", &self.loop_list.len())
            .field("loop_lookup", &self.loop_lookup)
            .field("num_pairs", &self.pair_list.len())
            .field("num_add_neighbors", &self.loop_neighbors.values().map(|v| v.len()).sum::<usize>())
            .field("num_del_neighbors", &self.pair_neighbors.len())
            .finish()
    }
}

impl<'a, M: EnergyModel> fmt::Display for LoopStructure<'a, M> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Convert sequence to string
        let mut dbr = vec!['.'; self.sequence.len()];
        for (i, j) in &self.pair_list {
            dbr[*i] = '(';
            dbr[*j] = ')';
        }
        let dbr_str: String = dbr.into_iter().collect();
        write!(f, "{}", dbr_str)
    }
}

pub fn allocate_index(indices: &mut IntSet<usize>, fallback: usize) -> usize {
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
 
    pub fn get_del_neighbors(&self) -> Vec<(usize, usize, i32)> {
        self.pair_neighbors
            .iter()
            .map(|(&i, &delta_e)| (i, self.pair_list[&i], delta_e))
            .collect()
    }

    pub fn energy(&self) -> i32 {
        self.loop_list
            .values()
            .map(|(_, e)| *e)
            .sum()
    }

    pub fn update_pair_neighbors(&mut self, pairs: &[(usize, usize)]
    ) -> Vec<(usize, usize, i32)> 
    {
        let mut change = vec![];
        for (i, j) in pairs {
            let (outer, outer_energy) = self.loop_list.get(self.loop_lookup.get(i).unwrap()).unwrap();
            let (inner, inner_energy) = self.loop_list.get(self.loop_lookup.get(j).unwrap()).unwrap();
            let combo = &outer.join_loop(&inner);
            let combo_energy = self.model.energy_of_loop(self.sequence, &combo);
            // How does the free energy change if the move is applied.
            let delta = combo_energy - (outer_energy + inner_energy);
            assert_eq!(self.pair_list[i], *j);
            self.pair_neighbors.insert(*i, delta);
            change.push((*i, *j, delta));
        }
        change
    }
  
    pub fn loop_neighbors(&self) -> &IntMap<usize, Vec<(usize, usize, i32)>> {
        &self.loop_neighbors
    }
 
    pub fn loop_lookup(&self) -> &IntMap<usize, usize> {
        &self.loop_lookup
    }

     /// Produce a new updated SecStruct.
    pub fn apply_del_move(&mut self, i: usize, j: usize) -> 
        ((usize, Vec<(usize, usize, i32)>), Vec<(usize, usize, i32)>
    ) {
        assert_eq!(&j, self.pair_list.get(&i).expect("Missing pair_list entry."));
        let &lli = self.loop_lookup.get(&i).expect("Missing loop_lookup entry for i.");
        let &llj = self.loop_lookup.get(&j).expect("Missing loop_lookup entry for j.");
        let (outer, o_en) = self.loop_list.get(&lli).expect("Missing loop_list entry lli.");
        let (inner, i_en) = self.loop_list.get(&llj).expect("Missing loop_list entry llj.");

        // Recompute or get from global cache.
        let combo = &outer.join_loop(&inner);
        let delta = self.pair_neighbors.get(&i).expect("Missing pair_neighbors entry.");
        let combo_energy = (o_en + i_en) + delta; // we are adding the change!
        //assert_eq!((o_en + i_en) + delta, self.model.energy_of_loop(self.sequence, &combo), "this should add up!");

        // Cleanup and re-allocate.
        self.pair_list.remove(&i); //todo
        self.pair_neighbors.remove(&i);//todo

        self.loop_list.remove(&lli);
        self.loop_list.remove(&llj);
        self.loop_neighbors.remove(&lli).expect("at least empty list.");
        self.loop_neighbors.remove(&llj).expect("at least empty list.");
        self.l_indices.insert(llj);

        // We re-use lli for the new loop.
        self.loop_list.insert(lli, (combo.clone(), combo_energy));
        let loop_neighbors = get_loop_neighbors(&combo, &combo_energy, self.sequence, self.model);
        self.loop_neighbors.insert(lli, loop_neighbors.clone());

        for k in &combo.inclusive_unpaired_indices(self.sequence.len()) {
            assert!(self.loop_lookup[k] == lli || self.loop_lookup[k] == llj);
            self.loop_lookup.insert(*k, lli);
        }

        let pair_changes = self.update_pair_neighbors(&combo.pairs());
        ((lli, loop_neighbors), pair_changes)
    }

    /// Produce a new updated SecStruct.
    pub fn apply_add_move(&mut self, i: usize, j: usize
    ) -> (
    (usize, Vec<(usize, usize, i32)>),
    (usize, Vec<(usize, usize, i32)>),
    Vec<(usize, usize, i32)>,
    ) {
        let model = self.model;
        let &lli = self.loop_lookup.get(&i).expect("Missing loop_lookup entry for i.");

        assert_eq!(lli, 
            *self.loop_lookup.get(&j).expect("Missing loop_lookup entry for j."),
            "Missing loop_lookup entry for j."
        );

        let (combo, c_energy) = {
            let (combo_ref, c_energy) = self.loop_list.get(&lli).expect("Missing loop_list entry lli.");
            (combo_ref.clone(), c_energy)
        };
        let (outer, inner) = combo.split_loop(i, j);
        //NOTE; could look delta up directly by searching loop_list.
        let outer_energy = model.energy_of_loop(self.sequence, &outer);
        let inner_energy = model.energy_of_loop(self.sequence, &inner);
        // How does the energy change if we apply the base-pair move.
        let delta = (outer_energy + inner_energy) - c_energy;

        self.loop_list.remove(&lli);
        self.loop_neighbors.remove(&lli).expect("there must have been at least an empty list");

        // First, let's reuse the lli we just freed.
        self.loop_list.insert(lli, (outer.clone(), outer_energy));
        let new_outer_add_neighbors = get_loop_neighbors(
            &outer, &outer_energy, self.sequence, model);
        let new_inner_add_neighbors = get_loop_neighbors(
            &inner, &inner_energy, self.sequence, model);
        self.loop_neighbors.insert(lli, new_outer_add_neighbors.clone());

        let llj = allocate_index(&mut self.l_indices, self.loop_list.len());
        self.loop_list.insert(llj, (inner.clone(), inner_energy));
        self.loop_neighbors.insert(llj, new_inner_add_neighbors.clone());

        self.pair_list.insert(i, j);

        for k in &outer.inclusive_unpaired_indices(self.sequence.len()) {
            self.loop_lookup.insert(*k, lli);
        }

        for k in &inner.inclusive_unpaired_indices(self.sequence.len()) {
            self.loop_lookup.insert(*k, llj);
        }

        // How does the energy change if we *revert* the base-pair move.
        self.pair_neighbors.insert(i, -delta);
        let mut pair_changes = self.update_pair_neighbors(&combo.pairs());
        pair_changes.extend(vec![(i, j, -delta)]);

        ((lli, new_outer_add_neighbors),
         (llj, new_inner_add_neighbors),
        pair_changes)
    }

}

impl<'a, T: LoopDecomposition, M: EnergyModel> TryFrom<(&'a [Base], &T, &'a M)> for LoopStructure<'a, M> {
    type Error = String;

    fn try_from((sequence, pairings, model): (&'a [Base], &T, &'a M)) -> Result<Self, Self::Error> {
        let mut loop_list = IntMap::default();
        let mut pair_list = IntMap::default();
        let mut loop_lookup = IntMap::default();
        let mut ll_idx = 0;
        //let mut pl_idx = 0;

        // Decomposing the structure into loops and initializing
        // loop_list, pair_list, and loop_lookup. 
        pairings.for_each_loop(|l| {
            // This it the point at which energy is evaluated and stored.
            let energy = model.energy_of_loop(sequence, l);
            loop_list.insert(ll_idx, (l.clone(), energy));
            if let Some((i, j)) = l.closing() {
                pair_list.insert(i, j); 
            }
            for k in &l.inclusive_unpaired_indices(sequence.len()) {
                loop_lookup.insert(*k, ll_idx);
            }
            ll_idx += 1;
        });

        // Now we want to get all neighbors where pairs can be added
        //
        let mut loop_neighbors = IntMap::default();
        for (nn_idx, (nn_loop, nn_energy)) in loop_list.iter() {
            let neighbors = get_loop_neighbors(&nn_loop, nn_energy, sequence, model);
            loop_neighbors.insert(*nn_idx, neighbors);
        }

        let mut pair_neighbors = IntMap::default();
        for (i, j) in pair_list.iter() {
            let (outer, outer_energy) = loop_list.get(loop_lookup.get(i).unwrap()).unwrap();
            let (inner, inner_energy) = loop_list.get(loop_lookup.get(j).unwrap()).unwrap();
            let combo = &outer.join_loop(&inner);
            let combo_energy = model.energy_of_loop(sequence, &combo);
            // How does the free energy change if the move is applied.
            let delta = combo_energy - (outer_energy + inner_energy);
            pair_neighbors.insert(*i, delta);
        }

        Ok(LoopStructure {
            sequence,
            model,
            loop_list,
            l_indices: IntSet::default(),
            pair_list,
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
        // How does the free energy change if the move is applied.
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


impl<'a, M: EnergyModel> From<&LoopStructure<'a, M>> for DotBracketVec {
    fn from(ls: &LoopStructure<'a, M>) -> Self {
        // Use the same logic as your Display impl, but avoid allocating a String unnecessarily
        let mut vec = vec![DotBracket::Unpaired; ls.sequence.len()];
        for (i, j) in &ls.pair_list {
            vec[*i] = DotBracket::Open;
            vec[*j] = DotBracket::Close;
        }
        DotBracketVec(vec)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use structure::PairTable;
    use energy::ViennaRNA;
    use energy::NucleotideVec;

    #[test]
    fn test_add_then_del_roundtrip() {
        let seq = NucleotideVec::from_lossy("GCCCCGGUCA");
        let structure = PairTable::try_from("...........").unwrap();
        let model = ViennaRNA::default();

        let mut ls = LoopStructure::try_from((&seq[..], &structure, &model)).unwrap();

        // Clone neighbor list so we don’t mutate while iterating
        let neighbors = ls.get_add_neighbors();

        for (i, j, de) in neighbors {
            let initial_energy = ls.energy();
            println!("({i} {j} {de}) at energy: {}", initial_energy);

            // add pair
            let _ = ls.apply_add_move(i, j);
            println!("{i} {j} {}", ls.energy());

            // delete the same pair
            let (p, q, rde) = ls.get_del_neighbors().first().cloned().unwrap();
            println!("({p} {q} {rde}) at energy: {}", ls.energy());
            assert_eq!((i, j), (p, q), "same pair gets deleted");
            assert_eq!(de, -rde, "inverse energy of reverse move");

            // delete pair 
            let _ = ls.apply_del_move(i, j);
            let roundtrip_energy = ls.energy();
            println!("{i} {j} {}", ls.energy());
            assert_eq!(roundtrip_energy, initial_energy, "roundtrip energy mismatch");
        }
    }

    #[test]
    fn test_add_then_del_bug() {
        let seq = NucleotideVec::from_lossy("GCCCCGGUCA");
        let structure = PairTable::try_from("((....).).").unwrap();
        let model = ViennaRNA::default();

        let ls = LoopStructure::try_from((&seq[..], &structure, &model)).unwrap();
        let neighbors = ls.get_del_neighbors();
        println!("{:?}", neighbors);

        let structure = PairTable::try_from("..........").unwrap();
        let mut ls = LoopStructure::try_from((&seq[..], &structure, &model)).unwrap();
        let _ = ls.apply_add_move(0, 8);
        println!("{:?}", neighbors);
        let _ = ls.apply_add_move(1, 6);
        assert_eq!(neighbors, ls.get_del_neighbors());
    }

}

