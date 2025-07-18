use domainlevel::Acfp;
use domainlevel::DomainRegistry;
use domainlevel::SegmentSequence;

fn main() {
    let acfp = Acfp::try_from(". () (). .(.) .(..) .((.))").unwrap();
    assert!(acfp.is_valid());

    let mut registry = DomainRegistry::new();
    let segseq = SegmentSequence::design_from_acfp(&acfp, &mut registry).unwrap();
    println!("{}", &segseq.get_domain_sequence().iter()
        .map(|d| format!("{}", d)).collect::<Vec<_>>().join(" "));

    for db in &segseq.get_domain_transcription_complexes(&registry) {
        println!("{:<2} {}", db.len(), db);
    }
    for db in &segseq.get_segment_transciption_complexes(&registry) {
        println!("{:<2} {}", db.len(), db);
    }
    for db in &segseq.get_adibatic_acfp(&registry) {
        println!("{:<2} {}", db.len(), db);
    }

    assert!(!segseq.implements_acfp(acfp.path(), &registry));
}


