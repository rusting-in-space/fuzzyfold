use std::io::{stdin, Write, BufRead};
use std::path::PathBuf;
use log::{debug};
use env_logger::Builder;
use clap::{Parser, ArgAction};
use anyhow::Result;
use anyhow::anyhow;

use energy::ViennaRNA;
use energy::EnergyModel;
use energy::NucleotideVec;
use structure::PairTable;
use structure::DotBracketVec;


/// ff-eval: Free energy evaluator
#[derive(Parser, Debug)]
#[command(name = "ff-eval")]
#[command(version, about = "Evaluate secondary structure free energies.")]
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

fn main() -> Result<()> {
    let cli = Cli::parse();
    init_logging(cli.verbose);
    let param_file = cli.model_parameters.unwrap_or_else(|| {
        PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/params/rna_turner2004.par"))
    });

    debug!("Using parameter file: {:?}", param_file);
    debug!("Temperature: {} °C", cli.temperature);
    let mut model = ViennaRNA::from_parameter_file(param_file)?;
    model.set_temperature(cli.temperature);


    let mut lines = stdin().lock().lines();

    // --- First line: sequence and maybe extra numbers ---
    let first = lines.next().ok_or_else(|| anyhow!("Missing first line"))??;

    let mut parts = first.split_whitespace();
    let token = parts.next().ok_or_else(|| anyhow!("Missing sequence"))?;
    let sequence = NucleotideVec::from_lossy(token);
    let _: Vec<f64> = parts.map(|s| s.parse().unwrap()).collect();
    println!("{}", sequence);

    // --- Remaining lines: dot-bracket + energy ---
    for line in lines {
        let line = line?;
        let mut parts = line.split_whitespace();

        let token = parts.next().ok_or_else(|| anyhow!("Missing structure"))?;
        let structure = DotBracketVec::try_from(token)?;
        let ref_en: f64 = parts.next().ok_or_else(|| anyhow!("Missing energy"))?.parse()?;

        let pairings = PairTable::try_from(&structure)?;
        let energy = model.energy_of_structure(&sequence, &pairings);

        let mark = if (ref_en * 100f64).round() as i32 != energy { "*" } else { "" };
        if true || mark == "*" {
            println!("{} {:6.2} {:6.2} {}", structure, energy as f64 / 100.0, ref_en, mark);
        }
    }

    Ok(())
}


