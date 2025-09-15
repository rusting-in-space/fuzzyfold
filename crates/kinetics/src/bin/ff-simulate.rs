use clap::Parser;
use anyhow;

use rand::rng;
use colored::*;
use structure::PairTable;
use energy::ViennaRNA;
use kinetics::LoopStructure;
use kinetics::LoopStructureSSA;
use kinetics::Metropolis;
use kinetics::commandline_utils::read_fasta_like_input;

/// SSA simulator for RNA folding
#[derive(Parser, Debug)]
#[command(name = "ff-simulate")]
#[command(version, about = "Stochastic Simulation Algorithm for RNA folding")]
struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    pub input: String,

    /// Base rate constant k0
    #[arg(long, default_value_t = 1.0)]
    k0: f64,

    /// Simulation start time
    #[arg(long, default_value_t = 0.0)]
    t0: f64,

    /// Simulation stop time
    #[arg(long, default_value_t = 1e7)]
    t8: f64,
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    // --- Build simulator ---
    let model = ViennaRNA::default();

    let (header, sequence, structure) = read_fasta_like_input(&cli.input)?;
    if let Some(h) = header {
        println!("{}", h.yellow())
    }

    let pairings = PairTable::try_from(&structure)?;
    let loops = LoopStructure::try_from((&sequence[..], &pairings, &model)).unwrap();

    let ratemodel = Metropolis::new(model.temperature(), cli.k0);
    let mut simulator = LoopStructureSSA::from((loops, &ratemodel));

    simulator.simulate(
        &mut rng(), 
        cli.t0, 
        cli.t8, 
        |t, ls| {
            println!("{:15.9} {} {:>6}", t, ls, ls.energy());
        }
    );
    Ok(())
}

