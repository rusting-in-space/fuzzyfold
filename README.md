# The fuzzyfold workspace

A rust-based open source collection of traditional nucleic acid folding algorithms.

ff-eval: (Re-)evaluate secondary structures, fuzzy ensemble free energies.
ff-simulate: Kinfold-like stochastic simulations.

# Student projects:
 - cotranscriptional folding of large RNAs (+ ALU analysis).
 - findpath/direct path sampling & barrier estimations.
 - smith-waterman / HMM-like read alignments.


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
This is a _very_ early stage coding project, which develops quickly. In other
words, you are welcome to use it for your research, but you may quickly
encounter its limitations. Furthermore, be prepared for frustrations from
drastic interface changes.

You may use GitHubs issues for suggestions, but we can also just schedule a
zoom call at this point. (The thought of an exploding issues section while 
there there are missing resources to adress them is unsetteling.)

# Developer Notes
Again, this is an early stage project, and "we" are still trying to figure out
a good system for this. Feel free to contact "us" and help with any part of
this project. Experience with other open-source code development projects would
be helpful, but is not a must.

# Roadmap to v1.0
 - EnergyModel interface (NearestNeighbor parameters, ViennaRNA, NUPACK compatibility).
 - Single-stranded nucleic acid folding kinetics (Kinfold-style).
 - Cotranscriptional folding and analysis.
 - Domain-level nucleic acid folding interface (dsdobjects, peppercorn).

# Notes / TODO
 - energy model
 - SSA
 - LoopStructure global cache interface. (after benchmarking)

