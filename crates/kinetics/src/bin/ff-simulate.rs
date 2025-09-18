use clap::Parser;
use energy::NucleotideVec;
use kinetics::plotting::plot_occupancy_over_time;
use kinetics::KineticModel;
use std::f64;
use std::path::PathBuf;
use clap::Args;
use std::sync::Arc;
use anyhow::Result;
use anyhow::bail;
use colored::*;

use rayon::prelude::*;

use kinetics::timeline::Timeline;
use kinetics::timeline::load_macrostates;
use rand::rng;
use structure::DotBracketVec;
use structure::PairTable;
use energy::ViennaRNA;
use energy::EnergyModel;
use kinetics::LoopStructure;
use kinetics::LoopStructureSSA;
use kinetics::Metropolis;
use kinetics::commandline_utils::read_fasta_like_input;
use energy::commandline_utils::EnergyModelArguments;

#[derive(Debug, Args)]
pub struct SimulationParameters {
    /// Base rate constant k0
    #[arg(long, default_value_t = 1e5)]
    pub k0: f64,

    /// Currently: The last time point of the linear scale.
    #[arg(long, default_value_t = 1e-5)]
    pub t_ext: f64,

    /// Simulation stop time (> t_ext).
    #[arg(long, default_value_t = 1.)]
    pub t_end: f64,

    /// Number of time points on the linear scale.
    #[arg(long, default_value_t = 1)]
    pub t_lin: usize,

    /// Number of time points on the logarithmic scale.
    #[arg(long, default_value_t = 20)]
    pub t_log: usize,
}

impl SimulationParameters {
    /// Validate that all parameters make sense.
    pub fn validate(&self) -> Result<()> {
        if self.t_end <= self.t_ext {
            bail!("t_end ({}) must be greater than t_ext ({})", self.t_end, self.t_ext);
        }
        if self.t_lin == 0 && self.t_log > 1 {
            bail!("t_lin must be > 0 if t_log > 1 (got t_lin={}, t_log={})", self.t_lin, self.t_log);
        }
        if self.k0 <= 0.0 {
            bail!("k0 must be positive (got {})", self.k0);
        }
        Ok(())
    }
}


#[derive(Debug, Args)]
pub struct IOParameters {
    /// Number of independent simulations
    #[arg(short, long, default_value_t = 1)]
    pub num_sims: usize,

    /// Optional macrostate definition files (each line = dot-bracket structure)
    #[arg(long, value_name = "FILE", num_args = 1.., required = false)]
    pub macrostates: Vec<PathBuf>,

    #[arg(long, default_value_t = false)]
    pub full_trajectory: bool,

    // These are the "debug modes"
    // full-trajectory (basic print, overrules all other output forms)
    // seeded-trajectory (provide seed or autogenerate)
    //    // mean-first-passage time (return actual mfpt of trajectory, not observed time.)
    // prints the mean fist passage time to any of the provided macrostates.
    // may assert that there is only one macrostate (with arbitraty many structures)
 
}


#[derive(Debug, Parser)]
#[command(name = "ff-simulate")]
#[command(version, about = "Stochastic Simulation Algorithm for RNA folding")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    #[command(flatten, next_help_heading = "I/O parameters")]
    pub io: IOParameters,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    pub simulation: SimulationParameters,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    pub energy: EnergyModelArguments,
}

fn get_output_times(
    t_ext: f64,
    t_end: f64,
    t_lin: usize,
    t_log: usize,
) -> Vec<f64> {
    let mut times = vec![0.0];

    // Linear segments: append `t_lin` evenly spaced points
    let start = *times.last().unwrap();
    let step = t_ext / t_lin as f64;
    for i in 1..=t_lin {
        times.push(start + i as f64 * step);
    }

    // Logarithmic tail
    let start = *times.last().unwrap();
    let log_start = start.ln();
    let log_end = t_end.ln();
    for i in 1..t_log {
        let frac = i as f64 / t_log as f64;
        let value = (log_start + frac * (log_end - log_start)).exp();
        times.push(value);
    }
    times.push(t_end);

    times
}

fn simu_callback_full_trajecory<E: EnergyModel, K: KineticModel>(
    sequence: &NucleotideVec, pairings: &PairTable, emodel: &E, rmodel: &K, t_end: f64
) {
    println!("{} {:>6} {:>12} {:>5}", sequence, "energy".green(), "timepoint".cyan(), "#ntp".cyan());
    let loops = LoopStructure::try_from((&sequence[..], pairings, emodel)).unwrap();
    let mut simulator = LoopStructureSSA::from((loops, rmodel));
    simulator.simulate(
        &mut rng(), 
        t_end, 
        |t, tinc, flux, ls| {
            println!("{} {:6.2} {:12.8} {} ", 
                ls,
                ls.energy(),
                t + tinc,
                flux
            );
        }
    );
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.simulation.validate()?;

    // --- Build simulator ---
    let model = ViennaRNA::default();
    let ratemodel = Metropolis::new(model.temperature(), cli.simulation.k0);

    let (header, sequence, structure) = read_fasta_like_input(&cli.input)?;
    let pairings = PairTable::try_from(&structure)?;
    if let Some(h) = header {
        println!("{}", h.yellow())
    }
    println!("{} {:>6} {:>12} {:>5}", sequence, "energy".green(), "timepoint".cyan(), "#ntp".cyan());

    if cli.io.full_trajectory {
        simu_callback_full_trajecory(&sequence, &pairings, &model, &ratemodel, cli.simulation.t_end);
        return Ok(());
    }

    let times = get_output_times(
        cli.simulation.t_ext, 
        cli.simulation.t_end,
        cli.simulation.t_lin, 
        cli.simulation.t_log);

    let registry = load_macrostates(
        &cli.io.macrostates,
        &sequence, 
        Some(&model));

    let mreg = Arc::new(registry);

    let timelines: Vec<Timeline> = (0..cli.io.num_sims)
        .into_par_iter()
        .map(|_| {
            let registry = Arc::clone(&mreg);
            let mut timeline = Timeline::new(&times, registry);

            let loops = LoopStructure::try_from((&sequence[..], &pairings, &model)).unwrap();
            let mut simulator = LoopStructureSSA::from((loops, &ratemodel));
            let mut t_idx = 0;
            simulator.simulate(&mut rng(), 
                cli.simulation.t_end,
                |t, tinc, _flux, ls| {
                    while t_idx < times.len() && t+tinc >= times[t_idx] {
                        let structure = DotBracketVec::from(ls);
                        timeline.assign_structure(t_idx, &structure);
                        t_idx += 1;
                    }
            });
            timeline
        })
        .collect();
    
    // Master timeline
    let mut master = Timeline::new(&times, Arc::clone(&mreg));

    // Merge results from all threads
    for timeline in timelines {
        master.merge(timeline);
    }

    println!("{}", master);
    let _ = plot_occupancy_over_time(&master, "myplot.svg", cli.simulation.t_ext, cli.simulation.t_end);

    Ok(())
}



