use std::fmt;
use std::fs;
use std::sync::Arc;
use anyhow::{Result, bail};
use nohash_hasher::IntMap;
use serde::{Serialize, Deserialize};
use structure::DotBracketVec; 

use crate::macrostates::MacrostateRegistry;

/// One time point with its ensemble of macrostates.
#[derive(Debug)]
pub struct Timepoint {
    /// Absolute time in seconds
    pub time: f64,
    /// Mapping from macrostate index → number of trajectories in this state
    pub ensemble: IntMap<usize, usize>,
    /// Total number of observations recorded at this timepoint
    pub counter: usize,
}

impl Timepoint {
    /// Create a new empty timepoint
    pub fn new(time: f64) -> Self {
        Self {
            time,
            ensemble: IntMap::default(),
            counter: 0,
        }
    }

    /// Add a count for the given macrostate index
    pub fn add(&mut self, macro_idx: usize) {
        *self.ensemble.entry(macro_idx).or_insert(0) += 1;
        self.counter += 1;
    }

    /// Get the count for a specific macrostate (or 0 if not present)
    pub fn count(&self, macro_idx: usize) -> usize {
        *self.ensemble.get(&macro_idx).unwrap_or(&0)
    }

    /// Return the occupancy (fraction of total) for a macrostate
    pub fn occupancy(&self, macro_idx: usize) -> f64 {
        if self.counter == 0 {
            0.0
        } else {
            self.count(macro_idx) as f64 / self.counter as f64
        }
    }

    /// Iterate over all macrostate counts
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.ensemble.iter().map(|(k, v)| (*k, *v))
    }

}

#[derive(Debug)]
pub struct Timeline {
    /// Registry of all macrostates (used to classify structures)
    pub registry: Arc<MacrostateRegistry>,

    /// One `Timepoint` per output time in the simulation
    pub points: Vec<Timepoint>,
}

impl Timeline {
    /// Build a new empty timeline for given times and an existing macrostate registry.
    pub fn new(times: &[f64], registry: Arc<MacrostateRegistry>) -> Self {
        let points = times.iter().map(|&t| Timepoint::new(t)).collect();
        Self { registry, points }
    }

    /// Classify a structure and add it to the timeline at the given time index.
    /// Returns `Some(macro_idx)` if classified, `None` if unclassified.
    pub fn assign_structure(&mut self, t_idx: usize, structure: &DotBracketVec) {
        let m_idx = self.registry.classify(structure);
        self.points[t_idx].add(m_idx);
    }

    /// Get a reference to a timepoint by index.
    pub fn point(&self, t_idx: usize) -> &Timepoint {
        &self.points[t_idx]
    }

    /// Iterate over all timepoints with their index.
    pub fn iter(&self) -> impl Iterator<Item = (usize, &Timepoint)> {
        self.points.iter().enumerate()
    }

    pub fn merge(&mut self, other: Timeline) {
        assert!(
            Arc::ptr_eq(&self.registry, &other.registry),
            "Cannot merge timelines with different registries"
        );
        assert_eq!(self.points.len(), other.points.len(),
        "Cannot merge timelines with different numbers of timepoints");

        for (self_tp, other_tp) in self.points.iter_mut().zip(other.points) {
            for (macro_idx, count) in other_tp.iter() {
                *self_tp.ensemble.entry(macro_idx).or_insert(0) += count;
            }
            self_tp.counter += other_tp.counter;
        }
    }
}

impl fmt::Display for Timeline {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // DRF header
        writeln!(f, "{:>13} {:>5} {:>12} {:>10} {:>25}", "time", "id", "occupancy", "energy", "macrostate")?;
        for tp in self.points.iter() {
            let time = tp.time;
            let total = tp.counter.max(1);

            // Collect ensemble into a vector so we can sort it
            let mut entries: Vec<_> = tp.iter().collect();

            // Sort by energy, None last
            entries.sort_by(|(a_idx, _), (b_idx, _)| {
                let e_a = self.registry.get(*a_idx).energy();
                let e_b = self.registry.get(*b_idx).energy(); 
                match (e_a, e_b) {
                    (Some(a), Some(b)) => a.partial_cmp(&b).unwrap_or(std::cmp::Ordering::Equal),
                    (Some(_), None) => std::cmp::Ordering::Less,
                    (None, Some(_)) => std::cmp::Ordering::Greater,
                    (None, None) => std::cmp::Ordering::Equal,
                }
            });

            // Sort ensemble by energy (you could make this configurable)
            for (m_idx, count) in entries {
                let occu = count as f64 / total as f64;

                let name = self.registry.get(m_idx).name();
                let energy = self.registry.get(m_idx).energy();

                writeln!(
                    f,
                    "{:13.9} {:5} {:12.8} {:>10} {:>25}",
                    time,
                    m_idx,
                    occu,
                    energy.map_or("N/A".to_string(), |e| format!("{:10.2}", e)),
                    name,
                )?;
            }
        }
        Ok(())
    }
}


#[derive(Serialize, Deserialize)]
pub struct SerializableTimeline {
    points: Vec<SerializableTimePoint>,
}

#[derive(Serialize, Deserialize)]
pub struct SerializableTimePoint {
    time: f64,
    ensemble: Vec<(String, usize)>, // (macrostate name, count)
    counter: usize,
}

impl Timeline {
    pub fn to_serializable(&self) -> SerializableTimeline {
        SerializableTimeline {
            points: self.points.iter().map(|tp| {
                let ensemble = tp.ensemble.iter()
                    .map(|(id, count)| {
                        let name = self.registry.get(*id).name().to_string();
                        (name, *count)
                    })
                    .collect();
                SerializableTimePoint {
                    time: tp.time,
                    ensemble,
                    counter: tp.counter,
                }
            }).collect()
        }
    }

    /// Load a timeline from a JSON file, checking against the provided registry
    pub fn from_file<P: AsRef<std::path::Path>>(
        path: P,
        times: &[f64],
        registry: Arc<MacrostateRegistry>,
    ) -> Result<Self> {
        let data = fs::read_to_string(path)?;
        let serial: SerializableTimeline = serde_json::from_str(&data)?;

        // Sanity check: number of timepoints must match
        if serial.points.len() != times.len() {
            bail!(
                "Mismatch: {} timepoints in file, but {} provided times",
                serial.points.len(),
                times.len()
            );
        }

        let mut timeline = Timeline::new(times, Arc::clone(&registry));

        for (tp, serial_tp) in timeline.points.iter_mut().zip(serial.points) {
            assert!(
                (tp.time - serial_tp.time).abs() < 1e-9,
                "Time mismatch at point: {} vs {}",
                tp.time,
                serial_tp.time
            );

            for (name, count) in serial_tp.ensemble {
                // Look up macrostate by name in registry
                if let Some((idx, _m)) = registry.iter().find(|(_, m)| m.name() == name) {
                    *tp.ensemble.entry(idx).or_insert(0) += count;
                    tp.counter += count;
                } else {
                    bail!("Macrostate '{}' not found in registry", name);
                }
            }
        }
        Ok(timeline)
    }
}


