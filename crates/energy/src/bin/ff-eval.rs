use env_logger::{Builder, Env};
use std::error::Error;
use std::io::{BufRead, stdin, Write};

use energy::basify;
use energy::ViennaRNA;
use energy::EnergyModel;
use structure::PairTable;

fn main() -> Result<(), Box<dyn Error>> {
    Builder::from_env(Env::default().default_filter_or("info"))
        .format(|buf, record| {
            writeln!(buf, "{}", record.args())
        })
    .init();
    let mut lines = stdin().lock().lines();
     

    // --- First line: sequence and maybe extra numbers ---
    let first = lines.next().ok_or("Missing first line")??;
    let mut parts = first.split_whitespace();
    let seq_str = parts.next().ok_or("Missing sequence")?;
    let seq = basify(seq_str);
    println!("{}", seq_str);

    // You may have trailing numbers after sequence:
    let _: Vec<f64> = parts.map(|s| s.parse().unwrap()).collect();

    // ToyModel params (tweak as you like)
    //let mut model = ViennaRNA::default();
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/params/rna_andronescu2007.par");
    let model = ViennaRNA::from_parameter_file(path).unwrap();

    // --- Remaining lines: dot-bracket + energy ---
    for line in lines {
        let line = line?;
        let mut parts = line.split_whitespace();

        let db = parts.next().ok_or("Missing structure")?;
        let energy: f64 = parts.next().ok_or("Missing energy")?.parse()?;

        let pt = PairTable::try_from(db)?;
        let myen = model.energy_of_structure(&seq, &pt);

        let mark = if (energy * 100f64).round() as i32 != myen { "*" } else { "" };
        if mark == "*" {
        println!("{}\n{} {} -> {} {}", seq_str, db, energy, myen as f64 / 100.0, mark);
        }
    }

    Ok(())
}

