# The fuzzyfold workspace

An open source collection of nucleic acid folding algorithms.

ff-eval: (Re-)evaluate secondary structures, fuzzy ensemble free energies.
ff-simulate: Kinfold-like stochastic simulations.

# Individual projects goals:
 - cotranscriptional folding of large RNAs (+ ALU analysis).
 - findpath/direct path sampling & barrier estimations.
 - smith-waterman / HMM-like read alignments.
 - ff-eval with new energy model & analysis.

Crates:
 - [master] fuzzychain.structure: A library of secondary structure representations.
 - [master] fuzzychain.energy: Models for secondary structure energy evaluation.
 - [master] fuzzychain.kinetics: Models for the kinetic folding of nucleic acids.
 - [development] fuzzychain.domainlevel: A library of secondary structure representations.

Branches: 
 - master: well-structured, well-documented, high-coverage, approved, production-ready code.
 - development: well-integrated, some documentation, some coverage, unapproved, experimental code.
 - feature-_name_: early code proposals.

# User Notes
This is a _very_ early stage, quickly developing, coding project. You are
welcome to use it for research, but be prepared for frustration from drastic
interface changes.

You may use GitHubs issues for suggestions, but we can also just schedule a
zoom call at this point. (The thought of an exploding issues section while 
there there are missing resources to adress them is unsetteling.)

# Developer Notes
Feel free to reach out and help with any part of this project. Experience with
other open-source code development projects would be helpful, but is not a
must.

# Roadmap to v1.0
 - EnergyModel interface (NearestNeighbor parameters, ViennaRNA, NUPACK compatibility).
 - Single-stranded nucleic acid folding kinetics (Kinfold-style).
 - Cotranscriptional folding and analysis.
 - Domain-level nucleic acid folding interface (dsdobjects, peppercorn).

# Work in progess:
 - energy model
 - SSA
 - LoopStructure global cache interface. (after benchmarking)


