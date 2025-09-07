use std::error::Error;
use std::io::{BufRead, stdin};

use energy::basify;
use energy::ViennaRNA;
use energy::EnergyModel;
use structure::PairTable;

fn main() -> Result<(), Box<dyn Error>> {
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
    let model = ViennaRNA::default();

    // --- Remaining lines: dot-bracket + energy ---
    for line in lines {
        let line = line?;
        let mut parts = line.split_whitespace();

        let db = parts.next().ok_or("Missing structure")?;
        let energy: f64 = parts.next().ok_or("Missing energy")?.parse()?;

        let pt = PairTable::try_from(db)?;
        let myen = model.energy_of_structure(&seq, &pt);
        println!("{} {} -> {}", db, energy, myen);
    }

    Ok(())
}

