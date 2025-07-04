
/// Checks whether two sequence symbols are complementary (e.g., "a" <-> "a*")
pub fn is_complement(a: &str, b: &str) -> bool {
    a == format!("{}*", b) || b == format!("{}*", a)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_complement() {
        assert!(is_complement("a", "a*"));
        assert!(is_complement("b*", "b"));
        assert!(!is_complement("a", "b*"));
        assert!(!is_complement("x", "x"));
    }
 
}
