use std::fs::File;
use std::io::{stdin, Write, BufRead, BufReader};
use std::path::PathBuf;
use colored::*; // brings in the `.red()`, `.blue()` etc.
use log::{info, debug};
use env_logger::Builder;
use clap::{Parser, ArgAction};
use anyhow;


use energy::basify;
use energy::ViennaRNA;
use energy::EnergyModel;
use structure::PairTable;


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

fn read_fasta_like(path: &str) -> anyhow::Result<(Option<String>, String, String)> {
    let reader: Box<dyn BufRead> = if path == "-" {
        Box::new(BufReader::new(stdin()))
    } else {
        Box::new(BufReader::new(File::open(path)?))
    };

    let mut header: Option<String> = None;
    let mut sequence: Option<String> = None;
    let mut structure: Option<String> = None;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            header = Some(line.to_string());
        } else if sequence.is_none() {
            sequence = Some(line.to_string());
            if let Some(c) = line.chars().find(|&c| matches!(c, '.' | '(' | ')')) {
                anyhow::bail!(
                    "Invalid character '{}' in sequence (no '.', '(' and ')' allowed)",
                    c 
                );
            }
        } else if structure.is_none() {
            structure = Some(line.to_string());
            if let Some(c) = line.chars().find(|&c| !matches!(c, '.' | '(' | ')')) {
                anyhow::bail!(
                    "Invalid character '{}' in structure (only '.', '(' and ')' allowed)",
                    c
                );
            }
            break;
        }
    }

    Ok((
        header,
        sequence.ok_or_else(|| anyhow::anyhow!("Missing sequence line"))?,
        structure.ok_or_else(|| anyhow::anyhow!("Missing structure line"))?,
    ))
}

fn ruler(len: usize) -> String {
    let mut s = String::new();
    for i in 0..=len {
        if i % 10 == 0 {
            s.push_str(&format!("{}", i / 10));
        } else if i % 10 == 5 {
            s.push(',');
        } else {
            s.push('.');
        }
    }
    s
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    init_logging(cli.verbose);

    let (header, sequence, structure) = read_fasta_like(&cli.input)?;
    let param_file = cli.model_parameters.unwrap_or_else(|| {
        PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/params/rna_turner2004.par"))
    });
    debug!("Using parameter file: {:?}", param_file);
    debug!("Temperature: {} °C", cli.temperature);

    if let Some(h) = header {
        println!("{}", h.yellow())
    }

    let mut model = ViennaRNA::from_parameter_file(param_file)?;
    model.set_temperature(cli.temperature);
    let seq = basify(&sequence);
    let pt = PairTable::try_from(&structure[..])?;
    let energy = model.energy_of_structure(&seq, &pt);
    info!("{}", ruler(sequence.len() - 1).magenta());
    println!("{}\n{} {}", sequence, structure, format!("{:>6.2}", energy as f64 / 100.0).green());
    info!("{}", ruler(sequence.len() - 1).magenta());

    //// You may have trailing numbers after sequence:
    //let _: Vec<f64> = parts.map(|s| s.parse().unwrap()).collect();

    //// ToyModel params (tweak as you like)
    ////let mut model = ViennaRNA::default();
    //let path = concat!(env!("CARGO_MANIFEST_DIR"), "/params/rna_andronescu2007.par");
    //let model = ViennaRNA::from_parameter_file(path).unwrap();

    //// --- Remaining lines: dot-bracket + energy ---
    //for line in lines {
    //    let line = line?;
    //    let mut parts = line.split_whitespace();

    //    let db = parts.next().ok_or("Missing structure")?;
    //    let energy: f64 = parts.next().ok_or("Missing energy")?.parse()?;

    //    let pt = PairTable::try_from(db)?;
    //    let myen = model.energy_of_structure(&seq, &pt);

    //    let mark = if (energy * 100f64).round() as i32 != myen { "*" } else { "" };
    //    if mark == "*" {
    //    println!("{}\n{} {} -> {} {}", seq_str, db, energy, myen as f64 / 100.0, mark);
    //    }
    //}

    Ok(())
}



//fn old_main() {
//    let args = Args::parse();
//
//    if atty::is(Stream::Stdin) {
//        println!("Enter dot-bracket strings (one per line). Provide empty line to finish:");
//        println!("....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8");
//    }
//    
//    if args.min_subtree_size < 4 {
//        panic!("Please use min-subtree-size >= 4!")
//    }
//
//    let stdin = io::stdin();
//    let mut dot_bracket_list = Vec::new();
//    for line in stdin.lock().lines() {
//        let s = match line {
//            Ok(s) => s,
//            Err(_) => break,
//        };
//        let s = s.trim();
//
//        if s.is_empty() {
//            break;
//        }
//
//        if !s.chars().all(|c| matches!(c, '(' | ')' | '.')) {
//            panic!("Invalid input: \"{}\" (only '.', '(' and ')' are allowed)", s);
//        }
//
//        dot_bracket_list.push(s.to_string());
//    }
// 
//    if let Some(c) = args.cutoff {
//        calculate_all_subtrees(&dot_bracket_list, c, 
//            args.min_subtree_size, 
//            args.min_subtree_size_difference,
//            args.min_dotbracket_size_difference);
//    } else {
//        calculate_all(&dot_bracket_list, None);
//    }
//}


