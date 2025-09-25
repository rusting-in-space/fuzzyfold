use std::fs::File;
use std::path::PathBuf;
use std::io::{BufRead, BufReader};
use ahash::AHashMap;

use structure::DotBracketVec;
use structure::PairTable;
use energy::NucleotideVec;
use energy::EnergyModel;

use crate::{K0, KB};

#[derive(Debug)]
pub enum Macrostate {
    /// Explicit registry of structures in the macrostate.
    Explicit(StructureRegistry),

    /// Placeholder for constraint-based macrostates (not yet implemented).
    Constraint(DotBracketVec),
}

impl Macrostate {
    pub fn name(&self) -> &str {
        match self {
            Macrostate::Explicit(registry) => {
                &registry.name()
            }
            Macrostate::Constraint(_) => {
                todo!("");
            }
        }
    }

    pub fn energy(&self) -> Option<f64> {
        match self {
            Macrostate::Explicit(registry) => {
                registry.energy()
            }
            Macrostate::Constraint(_) => {
                todo!("");
            }
        }
    }
}

#[derive(Debug)]
pub struct MacrostateRegistry {
    macrostates: Vec<Macrostate>,
}

impl MacrostateRegistry {
    // -> change to default??
    pub fn new() -> Self {
        let unassigned = Macrostate::Explicit(StructureRegistry {
            name: "Unassigned".to_string(),
            pool: Vec::new(),
            lookup: AHashMap::new(),
            energy: None,
        });

        Self {
            macrostates: vec![unassigned],
        }
    }

    pub fn from_files<M: EnergyModel>(
        files: &[PathBuf],
        sequence: &NucleotideVec,
        model: Option<&M>,
    ) -> Self {
        let mut registry = MacrostateRegistry::new(); // default.

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

            let file_seq = NucleotideVec::from_lossy(&seq_line.trim());
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

            // Precompute RT once if we have a model
            let mut reg = StructureRegistry::from_list(&name, &structures);

            if let Some(model) = model {
                let rt = KB * (K0 + model.temperature());
                reg.assign_energy(sequence, model, rt)
            }

            registry.insert(Macrostate::Explicit(reg));
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
            match ms {
                Macrostate::Explicit(reg) => {
                    if reg.contains(structure) {
                        matches.push(i);
                    }
                }
                Macrostate::Constraint(_pattern) => {
                    // Not implemented yet
                }
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

    /// Iterate over all macrostates
    pub fn iter(&self) -> impl Iterator<Item = (usize, &Macrostate)> {
        self.macrostates.iter().enumerate()
    }
}

#[derive(Debug, Default)]
pub struct StructureRegistry {
    name: String,
    pool: Vec<DotBracketVec>,
    lookup: AHashMap<DotBracketVec, usize>,
    energy: Option<f64>,
}

impl StructureRegistry {
    pub fn from_list(name: &str, structures: &[DotBracketVec]) -> Self {
        let mut pool = Vec::with_capacity(structures.len());
        let mut lookup = AHashMap::with_capacity(structures.len());

        for (i, dbv) in structures.iter().enumerate() {
            pool.push(dbv.clone());
            lookup.insert(dbv.clone(), i);
        }

        Self {
            name: name.to_owned(),
            pool,
            lookup,
            energy: None,
        }
    }

    pub fn assign_energy<M: EnergyModel>(&mut self, sequence: &NucleotideVec, model: &M, rt: f64,) {
        let mut q_sum = 0.0;
        for dbv in self.pool.iter() {
            let pt = PairTable::try_from(dbv)
                .expect("Invalid dot-bracket for energy evaluation");
            let en = model.energy_of_structure(sequence, &pt) as f64 / 100.0;
            q_sum += (-en / rt).exp();
        }
        self.energy = Some(-rt * q_sum.ln());
    }


    /// Return the name of this registry (macrostate label)
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Return the name of this registry (macrostate label)
    pub fn energy(&self) -> Option<f64> {
        self.energy
    }

    /// Check whether a structure is contained in this registry.
    pub fn contains(&self, structure: &DotBracketVec) -> bool {
        self.lookup.contains_key(structure)
    }

    /// Return the number of unique structures in this registry.
    pub fn len(&self) -> usize {
        self.pool.len()
    }

    /// Return an iterator over all stored structures.
    pub fn iter(&self) -> impl Iterator<Item = &DotBracketVec> {
        self.pool.iter()
    }
}

