use std::path::PathBuf;
use clap::Args;

/// Free energy evaluation parameters.
#[derive(Debug, Args)]
pub struct EnergyModelArguments {
    /// Temperature in Celsius
    #[arg(short, long, default_value = "37.0")]
    pub temperature: f64,

    /// Parameter file (e.g. rna_turner2004.par)
    #[arg(short, long, value_name = "FILE")]
    pub model_parameters: Option<PathBuf>,
}

impl EnergyModelArguments {
    /// Return the parameter file path, falling back to crate-relative default.
    pub fn param_file(&self) -> PathBuf {
        self.model_parameters.clone().unwrap_or_else(|| {
            PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/params/rna_turner2004.par"))
        })
    }
}


