use core::f64;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::PathBuf;
use ahash::AHashMap;

use ahash::AHashSet;
use ff_structure::DotBracketVec;
use ff_structure::PairTable;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;
use ff_energy::Base;

use crate::reaction::ApplyMove;
use crate::LoopStructure;
use crate::RateModel;
use crate::{K0, KB};

// A trait macrostate would not be able to store information.
// It makes more sense to support different *From* traits.
#[derive(Debug)]
pub struct Macrostate {
    name: String,
    ensemble: AHashMap<DotBracketVec, (i32, f64)>,
    ensemble_energy: Option<f64>,
}

impl Macrostate {
    pub fn from_list<'a, E: EnergyModel>(
        name: &str, 
        sequence: &'a [Base],
        structures: &[DotBracketVec], 
        energy_model: &'a E, 
    ) -> Self {
        let mut ensemble = AHashMap::with_capacity(structures.len());
        let rt = KB * (K0 + energy_model.temperature());

        let mut q_sum = 0.0;
        for dbv in structures {
            let pt = PairTable::try_from(dbv)
                .expect("Invalid dot-bracket for energy evaluation");
            let en = energy_model.energy_of_structure(sequence, &pt);
            //TODO: overflow corrections?
            let q = (-en as f64 / 100.0 / rt).exp();
            ensemble.insert(dbv.clone(), (en, q));
            q_sum += q;
        }
        // Normalize probabilities
        for (_dbv, (_en, prob)) in ensemble.iter_mut() {
            *prob /= q_sum;
        }

        Self {
            name: name.to_owned(),
            ensemble,
            ensemble_energy: Some(-rt * q_sum.ln()),
        }
    }

    pub fn from_macrostate_neighbors<'a, E: EnergyModel, R: RateModel>(
        macrostate: &Macrostate,
        sequence: &'a [Base],
        energy_model: &'a E, 
        rate_model: &'a R, 
    ) -> Self {

        let origin: AHashSet<DotBracketVec> = macrostate.ensemble.keys().cloned().collect();
        let mut visited = AHashSet::default();
        let mut ensemble = AHashMap::default();

        for dbr in &origin {
            find_neighbors::<E, R>(
                dbr,
                None,
                sequence,
                energy_model,
                rate_model,
                &origin,
                &mut visited,
                &mut ensemble,
            );
        }

        let mut k_sum = 0.0;
        for (dbv, (en, k_in)) in &ensemble {
            let pt = PairTable::try_from(dbv)
                .expect("Invalid dot-bracket for energy evaluation");
            debug_assert_eq!(*en, energy_model.energy_of_structure(sequence, &pt));
            //TODO: overflow corrections?
            k_sum += k_in;
        }
        // Normalize probabilities
        for (_dbv, (_en, k_in)) in ensemble.iter_mut() {
            *k_in /= k_sum;
        }

        Self {
            name: macrostate.name.to_owned(),
            ensemble,
            ensemble_energy: None,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns all secondary structures, their corresponding energy and probablility
    pub fn ensemble(&self) -> &AHashMap<DotBracketVec, (i32, f64)> {
        &self.ensemble
    }

    pub fn ensemble_energy(&self) -> Option<f64> {
        self.ensemble_energy
    }
    
    pub fn contains(&self, structure: &DotBracketVec) -> bool {
        self.ensemble.contains_key(structure)
    }

    //fn get_random_structure(&self) -> DotBracketVec;
    //fn get_probability(&self, structure: &DotBracketVec) -> f64;
}

fn find_neighbors<'a, E: EnergyModel, R: RateModel>(
    dbr: &DotBracketVec,
    lss_opt: Option<&LoopStructure<'a, E>>,
    sequence: &'a [Base],
    energy_model: &'a E,
    rate_model: &'a R,
    origin: &AHashSet<DotBracketVec>,
    visited: &mut AHashSet<DotBracketVec>,
    neighbors: &mut AHashMap<DotBracketVec, (i32, f64)>,
) {
    // Stop if already processed
    if !visited.insert(dbr.clone()) {
        return;
    }

    // So the outer loop does not produce 
    // all loop structures unnecessarily.
    let mut lss = match lss_opt {
        Some(existing) => existing.clone(),
        None => {
            let pairings = PairTable::try_from(dbr).unwrap();
            LoopStructure::try_from((sequence, &pairings, energy_model)).unwrap()
        }
    };

    let mut mdbr = dbr.clone();
    for (bp_move, delta) in lss.all_moves() {
        // Let's look up whether the move is 
        // worth applying based on DotBracketVec.
        mdbr.apply_move(bp_move);

        if origin.contains(&mdbr) {
            lss.apply_move(bp_move);
            find_neighbors(&mdbr, Some(&lss), sequence, energy_model, rate_model, origin, visited, neighbors);
            lss.undo_move(bp_move);
        } else {
            // First step outside the macrostate
            let rate = rate_model.rate(delta);
            neighbors
                .entry(mdbr.clone())
                .and_modify(|(e, k)| {
                    debug_assert!(*e == lss.energy() + delta);
                    *k += rate;
                })
            .or_insert((lss.energy() + delta, rate));
        }
        mdbr.undo_move(bp_move);
    }
}

/// A registy to collect macrostate definitions.
/// MacrostateRegistry::from((seq, emodel))
pub struct MacrostateRegistry<'a, E: EnergyModel> {
    sequence: &'a NucleotideVec,
    energy_model: &'a E,
    // by convention: macrostates[0] = unassigned.
    macrostates: Vec<Macrostate>,
    _neighbors: Option<Vec<Macrostate>>,
}

impl<'a, E: EnergyModel> From<(&'a NucleotideVec, &'a E)> for MacrostateRegistry<'a, E> {
    fn from((sequence, energy_model): (&'a NucleotideVec, &'a E)) -> Self {
        let macrostates = vec![Macrostate {
            name: "Unassigned".to_string(),
            ensemble: AHashMap::new(),
            ensemble_energy: None,
        }];

        MacrostateRegistry {
            sequence,
            energy_model,
            macrostates,
            _neighbors: None,
        }
    }
}

impl<'a, E: EnergyModel> MacrostateRegistry<'a, E> {

    //TODO: preliminary interface: directly from input files.
    pub fn insert_files(&mut self, files: &[PathBuf]) {
        for file in files {
            let fh = File::open(file)
                .unwrap_or_else(|_| panic!("Failed to open macrostate file {:?}", file));
            let mut lines = BufReader::new(fh).lines();

            // --- Parse header line (must exist)
            let header_line = lines
                .next()
                .expect("Macrostate file is missing a header line")
                .expect("Failed to read header line");

            let name = header_line
                .strip_prefix('>')
                .unwrap_or_else(|| panic!("First line must start with '>' in {:?}", file))
                .trim()
                .to_string();

            // --- Parse sequence line (must exist)
            let seq_line = lines
                .next()
                .expect("Macrostate file is missing a sequence line")
                .expect("Failed to read sequence line");

            let file_seq = NucleotideVec::from_lossy(seq_line.trim());
            assert_eq!(
                file_seq, *self.sequence,
                "Sequence in macrostate file {:?} does not match provided input sequence",
                file
            );

            // --- Parse structures (must be at least one)
            let mut structures = Vec::new();
            for (lineno, line) in lines.enumerate() {
                let line = line.expect("Failed to read line");
                let line = line.trim();
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                match DotBracketVec::try_from(line) {
                    Ok(dbv) => structures.push(dbv),
                    Err(e) => panic!(
                        "Error parsing dot-bracket at line {} in {:?}: {:?}",
                        lineno + 3, // +3 to account for header+sequence lines
                        file,
                        e
                    ),
                }
            }

            assert!(
                !structures.is_empty(),
                "Macrostate file {:?} does not contain any structures", file
            );

            self.macrostates.push(Macrostate::from_list(
                    &name, self.sequence, &structures, self.energy_model
            ));
        }
    }

    /// Try to classify a structure:
    /// - Returns Some(index) if exactly one macrostate contains the structure
    /// - Returns None if no macrostate matches
    /// - Panics if more than one macrostate matches (unimplemented)
    pub fn classify(&self, structure: &DotBracketVec) -> usize {
        let mut matches = Vec::new();

        for (i, ms) in self.macrostates.iter().enumerate() {
            if ms.contains(structure) {
                matches.push(i);
            }
        }

        match matches.len() {
            0 => 0,
            1 => matches[0],
            _ => {
                // For now: raise a panic, since overlapping macrostates are ambiguous
                panic!("Structure {:?} belongs to multiple macrostates — not implemented", structure);
            }
        }
    }

    /// Get a macrostate by index
    pub fn get(&self, idx: usize) -> &Macrostate {
        &self.macrostates[idx]
    }

    /// Number of macrostates, including the catch-all unassigned macrostate.
    pub fn len(&self) -> usize {
        self.macrostates.len()
    }

    //TODO: Useless: there is always one.
    pub fn is_empty(&self) -> bool {
        self.macrostates.is_empty()
    }

    /// Iterate over all macrostates
    pub fn iter(&self) -> impl Iterator<Item = (usize, &Macrostate)> {
        self.macrostates.iter().enumerate()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ff_energy::ViennaRNA;
    use crate::Metropolis;

    #[test]
    fn test_macrostate_init() {
        /*        
        >lmin=lm3_bh=3.0
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....)))).((((........))))...............
        .((((....)))).((((.(....).))))...............
        .((((....))))..(((........)))................
        .((((....)))).((((.(.....)))))...............
        .(((......))).((((........))))...............
        ..(((....)))..((((........))))...............
        .(((......)))..(((........)))................
        .(((.(...)))).((((........))))...............
        */        

        let seq = NucleotideVec::from_lossy("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC");
        let energy_model = ViennaRNA::default();
        let rate_model = Metropolis::new(energy_model.temperature(), 1.0);

        let s1 = DotBracketVec::try_from(".((((....)))).((((........))))...............").unwrap();
        let s2 = DotBracketVec::try_from(".((((....)))).((((.(....).))))...............").unwrap();
        let s3 = DotBracketVec::try_from(".((((....))))..(((........)))................").unwrap();
        let s4 = DotBracketVec::try_from(".((((....)))).((((.(.....)))))...............").unwrap();
        let s5 = DotBracketVec::try_from(".(((......))).((((........))))...............").unwrap();
        let s6 = DotBracketVec::try_from("..(((....)))..((((........))))...............").unwrap();
        let s7 = DotBracketVec::try_from(".(((......)))..(((........)))................").unwrap();
        let s8 = DotBracketVec::try_from(".(((.(...)))).((((........))))...............").unwrap();

        let macrostate = Macrostate::from_list(
            "lmin=lm3_bh=3.0",
            &seq, 
            &[s1, s2, s3, s4, s5, s6, s7, s8], 
            &energy_model
        );
        println!("Macrostate '{}':", macrostate.name());
        println!("  Ensemble size: {}", macrostate.ensemble().len());

        let ensemble = macrostate.ensemble().clone();
        let mut ensemble: Vec<_> = ensemble.iter().collect();
        ensemble.sort_by_key(|(_, (energy, _))| *energy);
        for (dbr, (energy, prob)) in ensemble.iter() {
            println!("  {} -> E(s) = {energy}, P(s) = {prob:.4}", dbr);
        }

        assert_eq!(macrostate.ensemble().len(), 8);

        let neighbors = Macrostate::from_macrostate_neighbors(
            &macrostate, 
            &seq, 
            &energy_model, 
            &rate_model
        );
        println!("Neighbors '{}':", neighbors.name());
        println!("  Ensemble size: {}", neighbors.ensemble().len());

        let ensemble = neighbors.ensemble().clone();
        let mut ensemble: Vec<_> = ensemble.iter().collect();
        ensemble.sort_by_key(|(_, (energy, _))| *energy);
        for (dbr, (energy, prob)) in ensemble.iter() {
            println!("  {} -> E(s) = {energy}, P(s) = {prob:.4}", dbr);
        }

        assert_eq!(neighbors.ensemble().len(), 345);
    }
}

