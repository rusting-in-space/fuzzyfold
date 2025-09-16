use clap::Parser;
use clap::Args;
use anyhow::Result;
use anyhow::bail;
use colored::*;

use rand::rng;
use structure::PairTable;
use energy::ViennaRNA;
use kinetics::LoopStructure;
use kinetics::LoopStructureSSA;
use kinetics::Metropolis;
use kinetics::commandline_utils::read_fasta_like_input;
use energy::commandline_utils::EnergyModelArguments;

#[derive(Debug, Args)]
pub struct SimulationParameters {
    /// Base rate constant k0
    #[arg(long, default_value_t = 1e5)]
    k0: f64,

    /// Currently: The last time point of the linear scale.
    #[arg(long, default_value_t = 1e-5)]
    t_ext: f64,

    /// Simulation stop time (> t_ext).
    #[arg(long, default_value_t = 60.)]
    t_end: f64,

    /// Number of time points on the linear scale.
    #[arg(long, default_value_t = 1)]
    t_lin: usize,

    /// Number of time points on the logarithmic scale.
    #[arg(long, default_value_t = 20)]
    t_log: usize,
}

#[derive(Debug, Parser)]
#[command(name = "ff-simulate")]
#[command(version, about = "Stochastic Simulation Algorithm for RNA folding")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    pub input: String,

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


fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.simulation.t_lin == 0 && cli.simulation.t_log > 1 {
        bail!("Invalid combination: t_lin = 0, t_log > 1");
    }

    // --- Build simulator ---
    let model = ViennaRNA::default();

    let (header, sequence, structure) = read_fasta_like_input(&cli.input)?;
    if let Some(h) = header {
        println!("{}", h.yellow())
    }

    let pairings = PairTable::try_from(&structure)?;
    let loops = LoopStructure::try_from((&sequence[..], &pairings, &model)).unwrap();

    let ratemodel = Metropolis::new(model.temperature(), cli.simulation.k0);
    let mut simulator = LoopStructureSSA::from((loops, &ratemodel));

    println!("{} {:>6} {:>12} {:>5}", sequence, "energy".green(), "timepoint".cyan(), "#ntp".cyan());
    let times = get_output_times(
        cli.simulation.t_ext, 
        cli.simulation.t_end,
        cli.simulation.t_lin, 
        cli.simulation.t_log);

    let mut t_idx = 0;
    simulator.simulate(
        &mut rng(), 
        cli.simulation.t_end, 
        |t, ls| {
            while t_idx < times.len() && t >= times[t_idx] {
                println!("{} {:6.2} {:12.8} {:5} ", 
                    ls,
                    ls.energy(),
                    times[t_idx],
                    t_idx,
                );
                t_idx += 1;
            }
        }
    );
    Ok(())
}

