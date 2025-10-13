mod error;
mod dotbracket;
mod pair_table;
mod multi_pair_table;
mod loop_table;

pub use error::*;
pub use dotbracket::*;
pub use pair_table::*;
pub use multi_pair_table::*;
pub use loop_table::*;


/// We use u16 (0 to 65k), which is plenty for indexing positions on a nucleic
/// acid. If you ever want to change this, beware that some dependencies use
/// intmaps for pairs by casting a pair (u16, u16) into one u32. In particular,
/// ff_kinetics may *assume* that NAIDX is set to u16.
pub type NAIDX = u16;


