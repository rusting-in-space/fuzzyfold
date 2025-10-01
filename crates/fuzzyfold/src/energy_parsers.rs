use std::path::PathBuf;
use ff_energy::ViennaRNA;
use clap::Args;
use log::debug;

/// Free energy evaluation parameters.
#[derive(Debug, Args)]
pub struct EnergyModelArguments {
    /// Temperature in Celsius
    #[arg(short, long, default_value = "37.0")]
    pub temperature: f64,

    /// Parameter file (defaults to rna_turner2004.par)
    #[arg(short, long, value_name = "FILE")]
    pub model_parameters: Option<PathBuf>,
}

impl EnergyModelArguments {
    pub fn build_model(&self) -> ViennaRNA {
        debug!("Using parameter file: {:?}", self.model_parameters);
        debug!("Temperature: {} °C", self.temperature);
        let mut model = if let Some(path) = &self.model_parameters {
            ViennaRNA::from_parameter_file(path)
                .expect("Failed to load parameter file")
        } else {
            ViennaRNA::default()
        };
        model.set_temperature(self.temperature);
        model
    }
}


