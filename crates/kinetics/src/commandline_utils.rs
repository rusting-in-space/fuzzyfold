use std::fs::File;
use std::io::stdin;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Cursor;
use anyhow::Result;
use anyhow::anyhow;
use structure::DotBracketVec;

use energy::NucleotideVec;

pub fn read_fasta_like_input(s: &str) -> Result<(Option<String>, NucleotideVec, DotBracketVec)> {
    if s == "-" {
        read_fasta_like_stdin()
    } else {
        read_fasta_like_file(s)
    } 
}

pub fn read_fasta_like_string(s: &str) -> Result<(Option<String>, NucleotideVec, DotBracketVec)> {
    read_fasta_like(Cursor::new(s))
}

pub fn read_fasta_like_file<P: AsRef<std::path::Path>>(path: P) -> Result<(Option<String>, NucleotideVec, DotBracketVec)> {
    let reader = BufReader::new(File::open(path)?);
    read_fasta_like(reader)
}

pub fn read_fasta_like_stdin() -> Result<(Option<String>, NucleotideVec, DotBracketVec)> {
    let reader = BufReader::new(stdin());
    read_fasta_like(reader)
}

pub fn read_fasta_like<R: BufRead>(reader: R) -> Result<
(Option<String>, NucleotideVec, DotBracketVec)> {
    let mut header: Option<String> = None;
    let mut sequence: Option<NucleotideVec> = None;
    let mut structure: Option<DotBracketVec> = None;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            if sequence.is_some() && structure.is_some() {
                break;
            } else {
                continue;
            }
        }

        if line.starts_with('>') {
            header = Some(line.to_string());
        } else if sequence.is_none() {
            let token = line.split_whitespace().next().unwrap();
            sequence = Some(NucleotideVec::try_from(token)?);
        } else if structure.is_none() {
            let token = line.split_whitespace().next().unwrap();
            structure = Some(DotBracketVec::try_from(token)?);
            break;
        }
    }
    
    let sequence = sequence.ok_or_else(|| anyhow!("Missing sequence line"))?;
    let structure = structure.unwrap_or_else(|| {
        DotBracketVec::try_from(".".repeat(sequence.len()).as_str())
            .expect("Failed to construct open chain structure")
    });

    if sequence.len() != structure.len() {
        return Err(anyhow!(
            "Sequence length ({}) and structure length ({}) do not match",
            sequence.len(),
            structure.len()
        ));
    }

    Ok((header, sequence, structure))
}

pub fn ruler(len: usize) -> String {
    let mut s = String::new();
    let mut c = 0;
    for i in 0..=len {
        if i % 10 == 0 {
            let t = format!("{}", i / 10);
            c = t.len() - 1;
            s.push_str(&t);
            continue;
        } else if c > 0 {
            c -= 1;
            continue;
        }
        if i % 10 == 5 {
            s.push(',');
        } else {
            s.push('.');
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ruler() {
        assert_eq!(ruler(0), "0");
        assert_eq!(ruler(1), "0.");
        assert_eq!(ruler(5), "0....,");
        assert_eq!(ruler(10), "0....,....1");
        assert_eq!(ruler(20), "0....,....1....,....2");
        assert_eq!(ruler(120).len(), 120 + 2);
    }

    #[test]
    fn test_read_fasta_like_basic() {
        let input = ">test\nACGAAAAAAA\n(((....)))\n";
        let (header, seq, dbv) = read_fasta_like_string(input).unwrap();

        assert_eq!(header, Some(">test".into()));
        assert_eq!(seq.to_string(), "ACGAAAAAAA");
        assert_eq!(dbv.to_string(), "(((....)))");

        let input = "ACGU\n....\n";
        let (header, seq, dbv) = read_fasta_like_string(input).unwrap();

        assert_eq!(header, None);
        assert_eq!(seq.to_string(), "ACGU");
        assert_eq!(dbv.to_string(), "....");
    }

    #[test]
    fn test_read_fasta_like_basic_invalid_input() {
        let input = ">bad\nACGU\n..X.\n";
        let res = read_fasta_like_string(input);
        assert!(res.is_err(), "Expected error for invalid structure");

        let input = ">missing\n(((....)))\n";
        let res = read_fasta_like_string(input);
        assert!(res.is_err());
    }

    #[test]
    fn test_read_fasta_like_length_mismatch() {
        let input = "\
>test_long
ACGUACGUACGUACGUACGU
((((....))))((((....))))
";
        let res = read_fasta_like_string(input);
        assert!(res.is_err(), "Expected error for length difference");

    }

}


