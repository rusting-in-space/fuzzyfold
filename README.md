# The fuzzyfold workspace

An open source collection of nucleic acid folding algorithms.

**Note**: This is a _very_ early stage, quickly developing, coding project. You
are welcome to use it for research, but be prepared for frustration from
drastic interface changes. You may use GitHubs issues for suggestions, but you
are also welcome to reach out directly at this point.

## Current software:
 - **ff-eval**: Free energy evaluation details for secondary structures.
 - **ff-tajectory**: Single stochastic nucleic acid folding trajectories.
 - **ff-timecourse**: Stochastic nucleic acid secondary structure ensemble simulations.

## Current crates:
 - fuzzychain.structure: Common secondary structure data structures.
 - fuzzychain.energy: Secondary structure free energy evaluation.
 - fuzzychain.kinetics: Stochastic folding kinetics for nucleic acids.

### Early stage crates:
 - fuzzychain.plotting: Utilities for visualization.
 - fuzzychain.domainlevel: Domain-level secondary structure kinetics.

## Developer notes:
Feel free to reach out and contribute to this project. Experience with other
open-source code development projects would be helpful, but is not a must.

For benchmarking of the stochastic simutation algorithm:

    ```cargo bench --bench stochastic_simulation```

### Branches: 
 - **master**: well-structured, well-documented, high-coverage, approved, production-ready code.
 - **development**: well-integrated, some documentation, some coverage, unapproved, experimental code.
 - **feature-_name_**: early code proposals.


