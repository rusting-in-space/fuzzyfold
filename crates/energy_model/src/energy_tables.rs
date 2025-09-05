
use std::fmt;
use std::fs::File;
use std::path::Path;
use std::io::{BufRead, BufReader};

use strum::EnumCount;
use strum_macros::EnumCount;
use rustc_hash::FxHashMap;

use crate::parameter_parsing::ParamFileSection;
use crate::parameter_parsing::SectionParser;

#[derive(Debug)]
pub enum ParamError {
    Io(std::io::Error),
    Parse(String),
    MissingValue(&'static str, usize),
    InvalidHairpinSize(usize),
}

impl std::error::Error for ParamError {}

impl From<std::io::Error> for ParamError {
    fn from(e: std::io::Error) -> Self {
        ParamError::Io(e)
    }
}

impl fmt::Display for ParamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParamError::Io(e) => write!(f, "I/O error: {}", e),
            ParamError::Parse(msg) => write!(f, "Parse error: {}", msg),
            ParamError::MissingValue(table, index) => {
                write!(f, "Missing value in parameter table '{}' at index {}", table, index)
            }
            ParamError::InvalidHairpinSize(n) => {
                write!(f, "Invalid hairpin size: {}", n)
            }
        }
    }
}

#[derive(Clone, Hash, Copy, Debug, Eq, PartialEq, EnumCount)]
pub enum Base { A, C, G, U, N }

impl TryFrom<char> for Base {
    type Error = ();
    fn try_from(c: char) -> Result<Self, ()> {
        Ok(match c.to_ascii_uppercase() {
            'A' => Base::A,
            'C' => Base::C,
            'G' => Base::G,
            'U' | 'T' => Base::U,
            'A'..='Z' => Base::N,
            _ => return Err(()),
        })
    }
}

pub fn basify(seq: &str) -> Vec<Base> {
    seq.chars()
        .map(Base::try_from)
        .collect::<Result<_, _>>()
        .unwrap_or_else(|_| panic!("Invalid character in sequence: {}", seq))
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, EnumCount)]
pub enum PairTypeRNA {
    AU,
    UA,
    CG,
    GC,
    GU,
    UG,
    NN,
}

const PAIR_LOOKUP: [[PairTypeRNA; Base::COUNT]; Base::COUNT] = {
    use Base::*;
    use PairTypeRNA::*;
    let mut table = [[NN; Base::COUNT]; Base::COUNT];
    table[A as usize][U as usize] = AU;
    table[U as usize][A as usize] = UA;
    table[C as usize][G as usize] = CG;
    table[G as usize][C as usize] = GC;
    table[G as usize][U as usize] = GU;
    table[U as usize][G as usize] = UG;
    table
};

pub fn pair_type(b5: Base, b3: Base) -> PairTypeRNA {
    PAIR_LOOKUP[b5 as usize][b3 as usize]
}

pub fn rev_pair_type(pt: &PairTypeRNA) -> PairTypeRNA {
    use PairTypeRNA::*;
    match pt {
        AU => UA,
        UA => AU,
        CG => GC,
        GC => CG,
        GU => UG,
        UG => GU,
        NN => NN,
    }
}

const P: usize = PairTypeRNA::COUNT;
const B: usize = Base::COUNT;

#[derive(Debug)]
pub struct EnergyTables {
    pub stack:            [[Option<i32>; P]; P],
    pub stack_enthalpies: [[Option<i32>; P]; P],
                                                                          
    pub mismatch_hairpin:            [[[Option<i32>; B]; B]; P],
    pub mismatch_hairpin_enthalpies: [[[Option<i32>; B]; B]; P],
    pub mismatch_interior:            [[[Option<i32>; B]; B]; P],
    pub mismatch_interior_enthalpies: [[[Option<i32>; B]; B]; P],
    pub mismatch_interior_1n:            [[[Option<i32>; B]; B]; P],
    pub mismatch_interior_1n_enthalpies: [[[Option<i32>; B]; B]; P],
    pub mismatch_interior_23:            [[[Option<i32>; B]; B]; P],
    pub mismatch_interior_23_enthalpies: [[[Option<i32>; B]; B]; P],
    pub mismatch_multi:            [[[Option<i32>; B]; B]; P],
    pub mismatch_multi_enthalpies: [[[Option<i32>; B]; B]; P],
    pub mismatch_exterior:            [[[Option<i32>; B]; B]; P],
    pub mismatch_exterior_enthalpies: [[[Option<i32>; B]; B]; P],
                                                                          
    pub dangle5:            [[Option<i32>; B]; P],
    pub dangle5_enthalpies: [[Option<i32>; B]; P],
    pub dangle3:            [[Option<i32>; B]; P],
    pub dangle3_enthalpies: [[Option<i32>; B]; P],
                                                                          
    pub int11:            Box<[[[[Option<i32>; B]; B]; P]; P]>,
    pub int11_enthalpies: Box<[[[[Option<i32>; B]; B]; P]; P]>,
    pub int21:            Box<[[[[[Option<i32>; B]; B]; B]; P]; P]>,
    pub int21_enthalpies: Box<[[[[[Option<i32>; B]; B]; B]; P]; P]>,
    pub int22:            Box<[[[[[[Option<i32>; B]; B]; B]; B]; P]; P]>,
    pub int22_enthalpies: Box<[[[[[[Option<i32>; B]; B]; B]; B]; P]; P]>,

    pub hairpin:            [Option<i32>; 31],
    pub hairpin_enthalpies: [Option<i32>; 31],
    pub bulge:            [Option<i32>; 31],
    pub bulge_enthalpies: [Option<i32>; 31],
    pub interior:            [Option<i32>; 31],
    pub interior_enthalpies: [Option<i32>; 31],

    pub ml_params: [Option<i32>; 6],
    pub ninio: [Option<i32>; 3],
    pub misc: [Option<i32>; 4],

    pub hairpin_sequences: FxHashMap<Vec<Base>, (i32, i32)>,
}

macro_rules! section_match {
    ($enum:expr, $line:expr, $tables:expr, $($struct:ident),+ $(,)?) => {
        match $enum {
            $(
                ParamFileSection::$struct(ref mut s) => s.parse_line($line, &mut $tables),
            )+
                _ => { },
        }
    };
}

impl EnergyTables {
    pub fn default() -> Self {
        EnergyTables {
            stack:            [[None; P]; P],
            stack_enthalpies: [[None; P]; P],

            mismatch_hairpin:            [[[None; B]; B]; P],
            mismatch_hairpin_enthalpies: [[[None; B]; B]; P],
            mismatch_interior:            [[[None; B]; B]; P],
            mismatch_interior_enthalpies: [[[None; B]; B]; P],
            mismatch_interior_1n:            [[[None; B]; B]; P],
            mismatch_interior_1n_enthalpies: [[[None; B]; B]; P],
            mismatch_interior_23:            [[[None; B]; B]; P],
            mismatch_interior_23_enthalpies: [[[None; B]; B]; P],
            mismatch_multi:            [[[None; B]; B]; P],
            mismatch_multi_enthalpies: [[[None; B]; B]; P],
            mismatch_exterior:            [[[None; B]; B]; P],
            mismatch_exterior_enthalpies: [[[None; B]; B]; P],
            dangle5:            [[None; B]; P],
            dangle5_enthalpies: [[None; B]; P],
            dangle3:            [[None; B]; P],
            dangle3_enthalpies: [[None; B]; P],

            int11:            Box::new([[[[None; B]; B]; P]; P]),
            int11_enthalpies: Box::new([[[[None; B]; B]; P]; P]),
            int21:            Box::new([[[[[None; B]; B]; B]; P]; P]),
            int21_enthalpies: Box::new([[[[[None; B]; B]; B]; P]; P]),
            int22:            Box::new([[[[[[None; B]; B]; B]; B]; P]; P]),
            int22_enthalpies: Box::new([[[[[[None; B]; B]; B]; B]; P]; P]),            

            hairpin: [None; 31],
            hairpin_enthalpies: [None; 31],
            bulge: [None; 31],
            bulge_enthalpies: [None; 31],
            interior: [None; 31],
            interior_enthalpies: [None; 31],
            ml_params: [None; 6],
            ninio: [None; 3],
            misc: [None; 4],

            hairpin_sequences: FxHashMap::default(),
        }
    }

    pub fn from_parameter_file<P: AsRef<Path>>(path: P) -> Result<Self, ParamError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        Self::from_reader(reader)
    }

    fn from_reader<R: BufRead>(reader: R) -> Result<Self, ParamError> {
        let mut tables = EnergyTables::default();
        let mut section = ParamFileSection::None;

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            if line.is_empty() || line.starts_with("/*"){
                continue;
            } 

            if let Some(rest) = line.strip_prefix("# ") {
                match ParamFileSection::try_from(rest.trim()) {
                    Ok(sec) => {
                        section = sec;
                    }
                    Err(_) => {
                        eprintln!("⚠️  Unknown parameter file section: {:?}", rest);
                        return Err(ParamError::Parse(rest.to_string()));
                    }
                }
                continue;
            }

            section_match!(section, line, tables, 
                Stack, 
                StackEnthalpies,
                MismatchHairpin,
                MismatchHairpinEnthalpies,
                MismatchInterior,
                MismatchInteriorEnthalpies,
                MismatchInterior1n,
                MismatchInterior1nEnthalpies,
                MismatchInterior23,
                MismatchInterior23Enthalpies,
                MismatchMulti,
                MismatchMultiEnthalpies,
                MismatchExterior,
                MismatchExteriorEnthalpies,
                Dangle5,
                Dangle5Enthalpies,
                Dangle3,
                Dangle3Enthalpies,
                Int11,
                Int11Enthalpies,
                Int21,
                Int21Enthalpies,
                Hairpin,
                HairpinEnthalpies,
                Bulge,
                BulgeEnthalpies,
                Interior,
                InteriorEnthalpies,
                MLParams,
                Ninio,
                Misc,
                HairpinSequences,
            );
        }
        Ok(tables) 
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_stack() {
        let dummy = r#"
# stack
/*  CG    GC    GU    UG    AU    UA    NN          */
  -240  -330  -210  -140  -210  -210  -140    /* CG */
  -330  -340  -250  -150  -220  -240  -150    /* GC */
  -210  -250   130   -50  -140  -130   130    /* GU */
  -140  -150   -50    30   -60  -100    30    /* UG */
  -210  -220  -140   -60  -110   -90   -60    /* AU */
  -210  -240  -130  -100   -90  -130   -90    /* UA */
  -140  -150   130    30   -60   -90   130    /* NN */
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();

        // Check one of the parsed entries
        assert_eq!(tables.stack[PairTypeRNA::CG as usize][PairTypeRNA::CG as usize], Some(-240));
        assert_eq!(tables.stack[PairTypeRNA::GC as usize][PairTypeRNA::CG as usize], Some(-330));
        assert_eq!(tables.stack[PairTypeRNA::GU as usize][PairTypeRNA::CG as usize], Some(-210));
    }

    #[test]
    fn test_parse_mismatch() {
        use Base::*;
        use PairTypeRNA::*;
        let dummy = r#"
# mismatch_hairpin                                                                                         
  -80  -100  -110  -100   -80    /* CG,E */   
 -140  -150  -150  -140  -150    /* CG,A */
  -80  -100  -110  -100   -80    /* CG,C */   
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();

        // Check one of the parsed entries
        assert_eq!(tables.mismatch_hairpin[CG as usize][N as usize][N as usize], Some(-80));
        assert_eq!(tables.mismatch_hairpin[CG as usize][N as usize][A as usize], Some(-100));
    }

    #[test]
    fn test_parse_dangle() {
        use Base::*;
        use PairTypeRNA::*;
        let dummy = r#"
# dangle5
/*   N     A     C     G     U          */
   -10   -50   -30   -20   -10    /* CG */
    -0   -20   -30    -0    -0    /* GC */
   -20   -30   -30   -40   -20    /* GU */
   -10   -30   -10   -20   -20    /* UG */
   -20   -30   -30   -40   -20    /* AU */
   -10   -30   -10   -20   -20    /* UA */
    -0   -20   -10    -0    -0    /* NN */
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();

        // Check one of the parsed entries
        assert_eq!(tables.dangle5[CG as usize][N as usize], Some(-10));
        assert_eq!(tables.dangle5[CG as usize][A as usize], Some(-50));
    }

    #[test]
    fn test_parse_int11() {
        use Base::*;
        use PairTypeRNA::*;
        let dummy = r#"
# int11
  90    90    50    50    50    /* CG,CG,N */
  90    90    50    50    50    /* CG,CG,A */
  50    50    50    50    50    /* CG,CG,C */
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();

        // Check one of the parsed entries
        assert_eq!(tables.int11[CG as usize][CG as usize][N as usize][N as usize], Some(90));
        assert_eq!(tables.int11[CG as usize][CG as usize][N as usize][A as usize], Some(90));
        assert_eq!(tables.int11[CG as usize][CG as usize][N as usize][C as usize], Some(50));
    }

    #[test]
    fn test_parse_int21() {
        use Base::*;
        use PairTypeRNA::*;
        let dummy = r#"
# int21
   230   230   230   230   230    /* CG,CG,N,N */
   230   230   230   230   230    /* CG,CG,N,A */
   230   230   230   230   230    /* CG,CG,N,C */
   230   230   230   230   230    /* CG,CG,N,G */
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();

        // Check one of the parsed entries
        assert_eq!(tables.int21[CG as usize][CG as usize][N as usize][N as usize][N as usize], Some(230));
        assert_eq!(tables.int21[CG as usize][CG as usize][N as usize][N as usize][A as usize], Some(230));
    }

    #[test]
    fn test_parse_loops() {
        let dummy = r#"
# hairpin
   INF   INF   INF   540   560   570   540   600   550   640
   650   660   670   680   690   690   700   710   710   720
   720   730   730   740   740   750   750   750   760   760
   770
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();

        // Check one of the parsed entries
        assert_eq!(tables.hairpin[0], None);
        assert_eq!(tables.hairpin[1], None);
        assert_eq!(tables.hairpin[2], None);
        assert_eq!(tables.hairpin[3], Some(540));
        assert_eq!(tables.hairpin[29], Some(760));
        assert_eq!(tables.hairpin[30], Some(770));
    }

    #[test]
    fn test_sequence_parsing() {
        let dummy = r#"
# Hexaloops
ACAGUACU     280   -1680
ACAGUGAU     360   -1140
ACAGUGCU     290   -1280
ACAGUGUU     180   -1540

# Tetraloops
CAACGG     550     690
CCAAGG     330   -1030
CCACGG     370    -330
CCCAGG     340    -890
CCGAGG     350    -660
CCGCGG     360    -750
CCUAGG     370    -350
CCUCGG     250   -1390
CUAAGG     360    -760
CUACGG     280   -1070
CUCAGG     370    -660
CUCCGG     270   -1290
CUGCGG     280   -1070
CUUAGG     350    -620
CUUCGG     370   -1530
CUUUGG     370    -680

# Triloops
CAACG     680    2370
GUUAC     690    1080
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();
        println!("{:?}", tables.hairpin_sequences);

        // Check one of the parsed entries
        assert_eq!(tables.hairpin_sequences[&basify("CCAAGG")], (330, -1030));
        assert_eq!(tables.hairpin_sequences[&basify("CAACG")], (680, 2370));
    }

}

