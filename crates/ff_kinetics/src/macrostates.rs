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

use crate::RateModel;
use crate::{K0, KB};

// A trait macrostate would not be able to store information.
// It makes more sense to support different from traits.
#[derive(Debug)]
pub struct Macrostate {
    name: String,
    ensemble: AHashMap<DotBracketVec, (i32, f64)>,
    ensemble_energy: f64,
    //NOTE: potentially using Once_cell for lazy initialization?
    neighbors: Option<AHashMap<DotBracketVec, (i32, f64)>>,
}

impl Macrostate {

    pub fn from_list<M: EnergyModel>(
        name: &str, 
        structures: &[DotBracketVec], 
        sequence: &NucleotideVec, 
        energy_model: &M, 
    ) -> Self {
        let mut ensemble = AHashMap::with_capacity(structures.len());
        let rt = KB * (K0 + energy_model.temperature());

        let mut q_sum = 0.0;
        for dbv in structures {
            let pt = PairTable::try_from(dbv)
                .expect("Invalid dot-bracket for energy evaluation");
            let en = energy_model.energy_of_structure(sequence, &pt);
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
            ensemble_energy: -rt * q_sum.ln(),
            neighbors: None,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns all secondary structures, their corresponding energy and probablility
    pub fn ensemble(&self) -> &AHashMap<DotBracketVec, (i32, f64)> {
        &self.ensemble
    }

    pub fn ensemble_energy(&self) -> f64 {
        self.ensemble_energy
    }
    
    pub fn contains(&self, structure: &DotBracketVec) -> bool {
        self.ensemble.contains_key(structure)
    }

    /// Returns all first steps out of the macrostate, 
    /// their corresponding energy and k_{ms,j}
    pub fn get_neighbors<R: RateModel>(&self, rate_model: &R
    ) -> &AHashMap<DotBracketVec, (i32, f64)> {
        //TODO for each structure enum neighbors,
        // if neighbor not in macrostate, add it to results,
        // if it does not exist, then update the energy and
        // assign rate, otherwise debug_assert energy and 
        // add the rate to the flux k_{ms, j}
        todo!("");
    }

    //fn get_random_structure(&self) -> DotBracketVec;
    //fn get_probability(&self, structure: &DotBracketVec) -> f64;
}

#[derive(Debug)]
pub struct MacrostateRegistry {
    macrostates: Vec<Macrostate>,
}

impl Default for MacrostateRegistry {
    fn default() -> Self {
        let unassigned = Macrostate {
            name: "Unassigned".to_string(),
            ensemble: AHashMap::default(),
            ensemble_energy: 0.0,
            neighbors: None,
        };

        Self {
            macrostates: vec![unassigned],
        }
    }
}

impl MacrostateRegistry {
    pub fn from_files<M: EnergyModel>(
        files: &[PathBuf],
        sequence: &NucleotideVec,
        model: &M,
    ) -> Self {
        let mut registry = MacrostateRegistry::default(); 

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
                file_seq, *sequence,
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
                "Macrostate file {:?} does not contain any structures",
                file
            );

            registry.insert(Macrostate::from_list(&name, &structures, sequence, model));
        }

        registry
    }

    pub fn insert(&mut self, macrostate: Macrostate) {
        self.macrostates.push(macrostate);
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

    /// Number of macrostates
    pub fn len(&self) -> usize {
        self.macrostates.len()
    }

    pub fn is_empty(&self) -> bool {
        self.macrostates.is_empty()
    }

    /// Iterate over all macrostates
    pub fn iter(&self) -> impl Iterator<Item = (usize, &Macrostate)> {
        self.macrostates.iter().enumerate()
    }
}

#[derive(Debug, Default)]
pub struct StructureRegistry {
    name: String,
    pool: AHashSet<DotBracketVec>,
}

impl StructureRegistry {
    pub fn from_list(name: &str, structures: &[DotBracketVec]) -> Self {
        let mut pool = AHashSet::with_capacity(structures.len());

        for (i, dbv) in structures.iter().enumerate() {
            pool.insert(dbv.clone());
        }

        Self {
            name: name.to_owned(),
            pool,
        }
    }

    /// Return the name of this registry (macrostate label)
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Return the number of unique structures in this registry.
    pub fn len(&self) -> usize {
        self.pool.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pool.is_empty()
    }

    /// Return an iterator over all stored structures.
    pub fn iter(&self) -> impl Iterator<Item = &DotBracketVec> {
        self.pool.iter()
    }
}

