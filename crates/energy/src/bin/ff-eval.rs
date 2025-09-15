use std::io::Write;
use log::{info, debug};
use colored::*;
use env_logger::Builder;
use clap::{Parser, Args, ArgAction};
use anyhow;

use energy::ViennaRNA;
use energy::EnergyModel;
use structure::PairTable;
use energy::commandline_utils::ruler;
use energy::commandline_utils::read_fasta_like_input;
use energy::commandline_utils::EnergyModelArguments;


#[derive(Debug, Args)]
pub struct EvalInput {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    pub input: String,

    /// Verbosity (-v = info, -vv = debug)
    #[arg(short, long, action = ArgAction::Count)]
    pub verbose: u8,
}


#[derive(Debug, Parser)]
#[command(name = "ff-eval")]
#[command(author, version, about)]
pub struct Cli {
    #[command(flatten)]
    pub eval: EvalInput,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    pub energy: EnergyModelArguments,
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
    init_logging(cli.eval.verbose);

    debug!("Using parameter file: {:?}", cli.energy.param_file());
    debug!("Temperature: {} °C", cli.energy.temperature);
    let mut model = ViennaRNA::from_parameter_file(cli.energy.param_file())?;
    model.set_temperature(cli.energy.temperature);

    //if atty::is(Stream::Stdin) {
    //    println!("Enter dot-bracket strings (one per line). Provide empty line to finish:");
    //    println!("....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8");
    //}

    let (header, sequence, structure) = read_fasta_like_input(&cli.eval.input)?;
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

