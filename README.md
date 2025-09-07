# The fuzzyfold workspace

An open source collection of nucleic acid folding algorithms.

**Note**: This is a _very_ early stage, quickly developing, coding project. You are
welcome to use it for research, but be prepared for frustration from drastic
interface changes.

You may use GitHubs issues for suggestions, but we can also just schedule a
zoom call at this point. (The thought of an exploding issues section while 
there there are missing resources to adress them is unsetteling.)

## Current software (in preparation):
 - **ff-eval**: (Re-)evaluate secondary structures, fuzzy ensemble free energies.
 - **ff-simulate**: Kinfold-like stochastic simulations.

### WS25 -- projects goals:
 - cotranscriptional folding of large RNAs (+ ALU analysis).
 - findpath/direct path sampling & barrier estimations.
 - smith-waterman / HMM-like read alignments.
 - ff-eval with new energy model & analysis.
 - local flooding algorithms & partition functions (barriers, ?)

## Current crates (in preparation):
 - fuzzychain.structure: A library of secondary structure data structures.
 - fuzzychain.energy: Secondary structure energy evaluation.
 - fuzzychain.kinetics: Kinetic folding of nucleic acids.

### Early stage crates:
 - fuzzychain.domainlevel: Domain-level secondary structure kinetics.
 - fuzzychain.plotting: Utilities for visualization.

# Developer notes:
Feel free to reach out and help with any part of this project. Experience with
other open-source code development projects would be helpful, but is not a
must.

## Branches: 
 - **master**: well-structured, well-documented, high-coverage, approved, production-ready code.
 - **development**: well-integrated, some documentation, some coverage, unapproved, experimental code.
 - **feature-_name_**: early code proposals.

## Roadmap to v1.0
 - EnergyModel interface (NearestNeighbor parameters, ViennaRNA, NUPACK compatibility).
 - Single-stranded nucleic acid folding kinetics (Kinfold-style).
 - Cotranscriptional folding and analysis.
 - Domain-level nucleic acid folding interface (dsdobjects, peppercorn).

### Work in progess:
 - energy model
 - SSA
 - LoopStructure global cache interface. (after benchmarking)


