use std::io::Write;
use std::path::PathBuf;
use log::{info, debug};
use colored::*;
use env_logger::Builder;
use clap::{Parser, ArgAction};
use anyhow;

use energy::ViennaRNA;
use energy::EnergyModel;
use structure::PairTable;
use energy::commandline_utils::ruler;
use energy::commandline_utils::read_fasta_like_input;

/// ff-eval: Free energy evaluator
#[derive(Parser, Debug)]
#[command(name = "ff-eval")]
#[command(version, about = "Evaluate secondary structure free energies from sequence/structure input")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    pub input: String,

    /// Temperature in Celsius
    #[arg(short, long, default_value = "37.0")]
    pub temperature: f64,

    /// Parameter file (e.g. rna_turner2004.par)
    #[arg(short, long, value_name = "FILE")]
    pub model_parameters: Option<PathBuf>,

    /// Verbosity (-v = info, -vv = debug)
    #[arg(short, long, action = ArgAction::Count)]
    pub verbose: u8,
}

fn init_logging(verbosity: u8) {
    let level = match verbosity {
        0 => "warn",
        1 => "info",
        _ => "debug",
    };

    Builder::from_env(env_logger::Env::default().default_filter_or(level))
        .format(|buf, record| {
            // no prefix, just the message
            writeln!(buf, "{}", record.args())
        })
        .init();
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    init_logging(cli.verbose);
    let param_file = cli.model_parameters.unwrap_or_else(|| {
        PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/params/rna_turner2004.par"))
    });

    debug!("Using parameter file: {:?}", param_file);
    debug!("Temperature: {} °C", cli.temperature);
    let mut model = ViennaRNA::from_parameter_file(param_file)?;
    model.set_temperature(cli.temperature);

    let (header, sequence, structure) = read_fasta_like_input(&cli.input)?;
    if let Some(h) = header {
        println!("{}", h.yellow())
    }

    let pairings = PairTable::try_from(&structure)?;
    let energy = model.energy_of_structure(&sequence, &pairings);

    info!("{}", ruler(sequence.len() - 1).magenta());
    println!("{}\n{} {}", sequence, structure, format!("{:>6.2}", energy as f64 / 100.0).green());
    info!("{}", ruler(sequence.len() - 1).magenta());

    Ok(())
}

