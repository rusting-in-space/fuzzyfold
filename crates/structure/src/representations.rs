
use std::fmt;

#[derive(Debug)]
pub enum StructureError {
    UnmatchedParenthesis,
    InvalidToken(String),
}

impl fmt::Display for StructureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            StructureError::UnmatchedParenthesis => write!(f, "Unmatched parenthesis!"),
            StructureError::InvalidToken(tok) => write!(f, "Invalid token encountered in kernel: '{}'", tok),
        }
    }
}


impl std::error::Error for StructureError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DotBracket {
    Unpaired, // '.'
    Open,     // '('
    Close,    // ')'
}

impl From<char> for DotBracket {
    fn from(c: char) -> Self {
        match c {
            '.' => DotBracket::Unpaired,
            '(' => DotBracket::Open,
            ')' => DotBracket::Close,
            _ => panic!("Invalid dot-bracket character: {}", c),
        }
    }
}

impl From<DotBracket> for char {
    fn from(db: DotBracket) -> Self {
        match db {
            DotBracket::Open => '(',
            DotBracket::Close => ')',
            DotBracket::Unpaired => '.',
        }
    }
}

pub struct DotBracketVec<'a>(pub &'a [DotBracket]);

impl fmt::Display for DotBracketVec<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for db in self.0 {
            write!(f, "{}", char::from(*db))?;
        }
        Ok(())
    }
}

pub fn make_pair_table(structure: &[DotBracket]) -> Vec<usize> {
    let mut table = vec![structure.len()];
    table.resize(structure.len() + 1, 0);

    let mut stack = Vec::new();
    for (i, &c) in structure.iter().enumerate() {
        match c {
            DotBracket::Open => stack.push(i + 1), // +1 because position 0 is length
            DotBracket::Close => {
                if let Some(j) = stack.pop() {
                    table[i + 1] = j;
                    table[j] = i + 1;
                }
            }
            _ => {}
        }
    }
    table
} 

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LoopInfo {
    Unpaired {l: usize },
    Paired {o: usize, i: usize},
}

pub fn make_loop_table(pair_table: &[usize]) -> Vec<LoopInfo> {
    let n = pair_table[0];
    let mut loop_table = vec![LoopInfo::Unpaired{ l: 0 }; n + 1];
    loop_table[0] = LoopInfo::Unpaired{ l: n };

    let mut loop_index = 0;
    let mut stack: Vec<(usize, usize)> = Vec::new();
    let mut mloop = 0;

    for i in 1..=n {
        let j = pair_table[i];
        if j == 0 {
            loop_table[i] = LoopInfo::Unpaired{ l: loop_index };
        } else if j > i {
            let outer_loop = loop_index;
            mloop += 1;
            loop_index = mloop;
            loop_table[i] = LoopInfo::Paired{ o: outer_loop, i: loop_index };
            stack.push((j, loop_index));
        } else if j < i {
            if let Some((_, inner_loop)) = stack.pop() {
                loop_index = stack.last().map(|&(_,l)|  l).unwrap_or(0);
                loop_table[i] = LoopInfo::Paired{ o: loop_index, i: inner_loop };
            } else {
                panic!("Unbalanced structure at position {i}");
            }
        }
    }
    loop_table
}

pub fn make_pair_table_from_pairs<I>(length: usize, pairs: I) -> Vec<usize> 
where I: IntoIterator<Item = (usize, usize)>,
{
    let mut table = vec![0; length + 1];
    table[0] = length;

    for (i, j) in pairs {
        // Convert 0-based to 1-based indexing
        let (i1, j1) = (i + 1, j + 1);
        table[i1] = j1;
        table[j1] = i1;
    }

    table
}

pub fn dotbracket_from_pairtable(pair_table: &[usize]) -> Vec<DotBracket> {
    let n = pair_table[0];
    let mut structure = vec![DotBracket::Unpaired; n];
    for i in 1..=n {
        let j = pair_table[i];
        if j > i {
            structure[i - 1] = DotBracket::Open;
            structure[j - 1] = DotBracket::Close;
        }
    }
    structure
}

/// Check if the substructure from i..j is well-formed:
/// - All pairings are internal to the interval
/// - `pair_table` uses 1-based indexing (ViennaRNA-style)
pub fn is_well_formed(pair_table: &[usize], i: usize, j: usize) -> bool {
    assert!(j <= pair_table.len()-1, "careful: expects 0-based input");
    for k in i..j {
        let l = pair_table[k + 1];
        if l != 0 {
            let l = l - 1; // convert to 0-based
            if l < i || l >= j {
                return false; // points outside
            }
        }
    }
    true
}

/// Checks whether two sequence symbols are complementary (e.g., "a" <-> "a*")
pub fn is_complement(a: &str, b: &str) -> bool {
    a == format!("{}*", b) || b == format!("{}*", a)
}

pub fn parse_kernel(kernel: &str) -> Result<(Vec<String>, Vec<DotBracket>), StructureError> {
    let tokens: Vec<&str> = kernel.split_whitespace().collect();
    let mut sequence: Vec<String> = Vec::new();
    let mut structure: Vec<DotBracket> = Vec::new();
    let mut stack: Vec<usize> = Vec::new();
    
    for token in tokens {
        if token == ")" {
            match stack.pop() {
                Some(open_idx) => {
                    let opener = &sequence[open_idx];
                    let inferred = if opener.ends_with('*') {
                        opener.trim_end_matches('*').to_string()
                    } else {
                        format!("{}*", opener)
                    };
                    sequence.push(inferred);
                    structure.push(DotBracket::Close);
                }
                None => return Err(StructureError::UnmatchedParenthesis),
            }
        } else if token.ends_with('(') {
            let seq = token.strip_suffix('(').unwrap();
            if seq.is_empty() {
                return Err(StructureError::InvalidToken(token.to_string()));
            }
            sequence.push(seq.to_string());
            structure.push(DotBracket::Open);
            stack.push(sequence.len() - 1);
        } else if token.ends_with(')') {
            return Err(StructureError::InvalidToken(token.to_string()));
        } else if token.contains('(') || token.contains(')') {
            return Err(StructureError::InvalidToken(token.to_string()));
        } else {
            sequence.push(token.to_string());
            structure.push(DotBracket::Unpaired)
        }
    }
    if !stack.is_empty() {
        return Err(StructureError::UnmatchedParenthesis);
    }
    Ok((sequence, structure))
}

pub fn get_kernel(sequence: &[String], structure: &[DotBracket]) -> String {
    assert_eq!(sequence.len(), structure.len(), "Mismatched sequence and structure lengths");
    let mut tokens = Vec::new();
    for (seq, struc) in sequence.iter().zip(structure.iter()) {
        let token = match struc {
            DotBracket::Unpaired => seq.clone(), // just the sequence
            DotBracket::Open => format!("{}(", seq), // explicit opener
            DotBracket::Close => ")".to_string(), // implicit sequence, inferred
        };
        tokens.push(token);
    }
    tokens.join(" ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_kernel_unpaired() {
        let (seq, struc) = parse_kernel("a b c").unwrap();
        assert_eq!(seq, vec!["a", "b", "c"]);
        assert_eq!(struc, vec![
            DotBracket::Unpaired,
            DotBracket::Unpaired,
            DotBracket::Unpaired,
        ]);
    }

    #[test]
    fn test_parse_kernel_with_pairing() {
        let (seq, struc) = parse_kernel("x( y )").unwrap();
        assert_eq!(seq, vec!["x", "y", "x*"]);
        assert_eq!(struc, vec![DotBracket::Open, DotBracket::Unpaired, DotBracket::Close]);
    }

    #[test]
    fn test_parse_kernel_nested() {
        let (seq, struc) = parse_kernel("a( b( c ) d )").unwrap();
        assert_eq!(seq, vec!["a", "b", "c", "b*", "d", "a*"]);
        assert_eq!(struc, vec![
            DotBracket::Open,
            DotBracket::Open,
            DotBracket::Unpaired,
            DotBracket::Close,
            DotBracket::Unpaired,
            DotBracket::Close,
        ]);
    }

    #[test]
    fn test_parse_kernel_errors() {
        assert!(matches!(parse_kernel("a("), Err(StructureError::UnmatchedParenthesis)));
        assert!(matches!(parse_kernel("a)"), Err(StructureError::InvalidToken(_))));
        assert!(matches!(parse_kernel("x((y)"), Err(StructureError::InvalidToken(_))));
    }

    #[test]
    fn test_make_pair_table_simple() {
        let struct_vec = vec![
            DotBracket::Open,
            DotBracket::Unpaired,
            DotBracket::Unpaired,
            DotBracket::Close,
        ];
        let pt = make_pair_table(&struct_vec);
        assert_eq!(pt, vec![4, 4, 0, 0, 1]);
    }

    #[test]
    fn test_get_kernel_roundtrip() {
        let (seq, struc) = parse_kernel("x( y z )").unwrap();
        let kernel = get_kernel(&seq, &struc);
        assert_eq!(kernel, "x( y z )");
    }

    #[test]
    fn test_is_complement() {
        assert!(is_complement("a", "a*"));
        assert!(is_complement("b*", "b"));
        assert!(!is_complement("a", "b*"));
        assert!(!is_complement("x", "x"));
    }

    #[test]
    fn test_well_formed_empty_interval() {
        let pt = vec![3, 0, 0, 0]; // 3 elements (1-based: positions 1..3)
        assert!(is_well_formed(&pt, 0, 0)); 
        assert!(is_well_formed(&pt, 0, 1)); 
        assert!(is_well_formed(&pt, 0, 2)); 
        assert!(is_well_formed(&pt, 0, 3)); 
        assert!(is_well_formed(&pt, 1, 3)); 
        assert!(is_well_formed(&pt, 2, 3)); 
        assert!(is_well_formed(&pt, 3, 3)); 
    }

    #[test]
    fn test_well_formed_pairings_within_interval() {
        let pt = vec![5, 0, 4, 0, 2, 0]; // 1 pairs with 3, 2 unpaired (1-based)
        assert!(is_well_formed(&pt, 0, 5)); // Full interval -- 0-based
        assert!(is_well_formed(&pt, 0, 4)); 
        assert!(is_well_formed(&pt, 1, 5));
        assert!(is_well_formed(&pt, 1, 4));
        assert!(is_well_formed(&pt, 1, 4));
        assert!(is_well_formed(&pt, 2, 3));
        assert!(!is_well_formed(&pt, 0, 3)); 
        assert!(!is_well_formed(&pt, 1, 3)); 
        assert!(!is_well_formed(&pt, 2, 4)); 
    }

    #[test]
    #[should_panic(expected = "careful: expects 0-based input")]
    fn test_well_formed_out_of_bounds_assert() {
        let pt = vec![2, 0, 0];
        is_well_formed(&pt, 0, 3); // j = pt.len(), should panic
    }

    fn format_loop_info(loop_info: &[LoopInfo]) -> String {
        let mut out = Vec::new();
        for info in loop_info.iter().skip(1) { // skip index 0 (length)
            let s = match info {
                LoopInfo::Unpaired { l } => format!("{}", l),
                LoopInfo::Paired { o, i } => format!("{}/{}", o, i),
            };
            out.push(s);
        }
        format!("[{}]", out.join(", "))
    }

    #[test]
    fn test_pair_table_to_loop_index_01() {
        use LoopInfo::*;
        let dbs = ".(((...)).((...))..(.(...)))";
        let li = vec![ Unpaired { l: 28 },
                       Unpaired { l: 0 }, 
                       Paired { o: 0, i: 1 }, Paired { o: 1, i: 2 }, Paired { o: 2, i: 3 }, 
                       Unpaired { l: 3 }, Unpaired { l: 3 }, Unpaired { l: 3 }, 
                       Paired { o: 2, i: 3 }, Paired { o: 1, i: 2 },
                       Unpaired { l: 1 }, 
                       Paired { o: 1, i: 4 }, Paired { o: 4, i: 5 }, 
                       Unpaired { l: 5 }, Unpaired { l: 5 }, Unpaired { l: 5 }, 
                       Paired { o: 4, i: 5 }, Paired { o: 1, i: 4 }, 
                       Unpaired { l: 1 }, Unpaired { l: 1 }, 
                       Paired { o: 1, i: 6 }, 
                       Unpaired { l: 6 }, 
                       Paired { o: 6, i: 7 }, 
                       Unpaired { l: 7 }, Unpaired { l: 7 }, Unpaired { l: 7 }, 
                       Paired { o: 6, i: 7 }, Paired { o: 1, i: 6 }, Paired { o: 0, i: 1 }
        ];
        let db: Vec<DotBracket> = dbs.chars()
            .map(|c| DotBracket::from(c)).collect();

        let pt = make_pair_table(&db);
        let re = make_loop_table(&pt);
        println!("{:?}", re);
        println!("{}", format_loop_info(&re));
        assert_eq!(re, li); 

    }

    #[test]
    fn test_pair_table_to_loop_index_02() {
        use LoopInfo::*;
        let dbs = ".(((...)(...).((.(...))).)).";
        let li = vec![ Unpaired { l: 28 }, 
                       Unpaired { l: 0 }, 
                       Paired { o: 0, i: 1 },
                       Paired { o: 1, i: 2 },
                       Paired { o: 2, i: 3 },
                       Unpaired { l: 3 },
                       Unpaired { l: 3 }, 
                       Unpaired { l: 3 }, 
                       Paired { o: 2, i: 3 }, 
                       Paired { o: 2, i: 4 }, 
                       Unpaired { l: 4 }, 
                       Unpaired { l: 4 }, 
                       Unpaired { l: 4 },
                       Paired { o: 2, i: 4 }, 
                       Unpaired { l: 2 }, 
                       Paired { o: 2, i: 5 }, 
                       Paired { o: 5, i: 6 }, 
                       Unpaired { l: 6 }, 
                       Paired { o: 6, i: 7 }, 
                       Unpaired { l: 7 }, 
                       Unpaired { l: 7 }, 
                       Unpaired { l: 7 }, 
                       Paired { o: 6, i: 7 }, 
                       Paired { o: 5, i: 6 }, 
                       Paired { o: 2, i: 5 }, 
                       Unpaired { l: 2 }, 
                       Paired { o: 1, i: 2 }, 
                       Paired { o: 0, i: 1 },
                       Unpaired { l: 0 }
        ];

        let db: Vec<DotBracket> = dbs.chars()
            .map(|c| DotBracket::from(c)).collect();
        let pt = make_pair_table(&db);
        let re = make_loop_table(&pt);
        println!("{:?}", re);
        println!("{}", format_loop_info(&re));
        assert_eq!(re, li); 
    }

    #[test]
    fn test_pair_table_to_loop_index_03() {
        use LoopInfo::*;
        let dbs = ".(((...)(...))).((((.(...))).)).";
        let li = vec![ Unpaired { l: 32 }, Unpaired { l: 0 },
                       Paired { o: 0, i: 1 }, Paired { o: 1, i: 2 }, Paired { o: 2, i: 3 }, 
                       Unpaired { l: 3 }, Unpaired { l: 3 }, Unpaired { l: 3 }, 
                       Paired { o: 2, i: 3 }, Paired { o: 2, i: 4 }, 
                       Unpaired { l: 4 }, Unpaired { l: 4 }, Unpaired { l: 4 }, 
                       Paired { o: 2, i: 4 }, Paired { o: 1, i: 2 }, Paired { o: 0, i: 1 }, 
                       Unpaired { l: 0 }, 
                       Paired { o: 0, i: 5 }, Paired { o: 5, i: 6 }, Paired { o: 6, i: 7 }, Paired { o: 7, i: 8 }, 
                       Unpaired { l: 8 }, 
                       Paired { o: 8, i: 9 },
                       Unpaired { l: 9 }, Unpaired { l: 9 }, Unpaired { l: 9 }, 
                       Paired { o: 8, i: 9 }, Paired { o: 7, i: 8 }, Paired { o: 6, i: 7 }, 
                       Unpaired { l: 6 }, Paired { o: 5, i: 6 }, Paired { o: 0, i: 5 }, 
                       Unpaired { l: 0 }
        ];
        
        let db: Vec<DotBracket> = dbs.chars()
            .map(|c| DotBracket::from(c)).collect();

        let pt = make_pair_table(&db);
        let re = make_loop_table(&pt);
        println!("{:?}", re);
        println!("{}", format_loop_info(&re));
        assert_eq!(re, li); 
    }

}

