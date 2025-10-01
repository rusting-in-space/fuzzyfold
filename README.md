# The fuzzyfold workspace

An open source collection of nucleic acid folding algorithms.

**Note**: This is a _very_ early stage, quickly developing, coding project. You
are welcome to use it for research, but be prepared for frustration from
drastic interface changes. You may use GitHubs issues for suggestions, but you
are also welcome to reach out directly at this point.

## Current fuzzyfold software (preliminary):
 - **ff-eval**: Free energy evaluation details for secondary structures.
 - **ff-tajectory**: Single stochastic nucleic acid folding trajectories.
 - **ff-timecourse**: Stochastic nucleic acid secondary structure ensemble simulations.
 - **ff-randseq**: Generate a random sequence. 

(Other software is work in progress and not published to crates.io)

## Current fuzzyfold crates (preliminary):
 - ff_stucture: Common secondary structure data structures.
 - ff_energy: Secondary structure free energy evaluation.
 - ff_kinetics: Stochastic folding kinetics for nucleic acids.

(Other crates are work in progress and not published to crates.io)

## Developer notes:
Feel free to reach out and contribute to this project. Experience with other
open-source code development projects would be helpful, but is not a must.

Goal of the fuzzyfold workspace is to provide both production-ready code for
scientic analysis, while also inviting code contributions from a larger
community. Thus, we welcome new developers to contribute to both new and
existing crates. However, the fuzzyfold crate itself should only re-export the
"highest tier" production-ready code and corresponding commandline interfaces.

For benchmarking of the stochastic simutation algorithm:

    ```cargo bench --workspace```

### Git branches (work in progress): 
 - **master**: well-structured, well-documented, high-coverage, publication-ready code.
 - **development**: well-integrated, some documentation, some coverage, experimental code.
 - **feature-_name_**: early code/crate proposals.


