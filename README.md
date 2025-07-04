# Vienna Nucleic Acid Suite

A collection of traditional ViennaRNA algorithms for nucleic acid thermodynamics and kinetics.

Current focus:
 - single-stranded nucleic acid folding (Kinfold)
 - cotranscriptional folding 
 - cotranscriptional folding visualization.
 - domain-level reaction network enumeration.


## TODO / Outline:

* crates/structure/error.rs
* crates/structure/dotbracket.rs
* crates/structure/pair_table.rs
* crates/structure/loop_index.rs
* crates/structure/pair_list.rs

* crates/domainlevel/complex.rs
* crates/domainlevel/utils.rs

next steps: 
- crates/domainlevel/rules/
- crates/domainlevel/reaction.rs
- crates/domainlevel/acfps.rs

consider:
- crates/structure/ring_list.rs
- crates/structure/tree.rs

- crates/reactions/network.rs

utils/parsing.rs
utils/cycledetection.rs
utils/nussinov.rs
acfp/
