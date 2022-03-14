# nk_walk.cpp 

 ##### Main routine to simulate adative walks on the Stuart Kauffman's NK-fitness landscape model, with point mutations (as usual) and inversions (as new).


### Compile nk_walk.cpp and run the simulation as 
 ./nk_walk  -n **N** -k **K** -a **A** -e **epi** -snk **nk_seed** -swlk **wlk_seed** -t **T** -m **M** -f  **genome.csv** 

where, 
- **N**: Genome length.
- **K**: Number of epistatic interactions (0 <= **K** < **N**).
- **A**: Alphabet size (default = 2).
- **epi**: Type of epistatic interactions ('**ADJ**' (default) or '**RND**').
- **nk_seed**: PRNG seed value for the landscape (default = -1).
- **wlk_seed**: PRNG seed value for the random walk (default = -1).
- **T**: Number of steps in the random walk (default = 10000).
- **M**: Mutation type (float) (0.0 <= **M** <= 1.0)," such that **M** = 0.0: point mutations only (default)," **M** = 1.0: inversions only."
- **genome.csv**: initial genome (if no file containing the genome is indicated, then the initial genome is randomly generated).  

 Help: run the command line ./nk_walk -help  to print out the required 
 input parameters. 

 #### Why reinvent the wheel? 
 This code uses the one developed by Wim Hordijk (in its version of 
 August 23, 2010 and which is available at http://www.cs.unibo.it/~fioretti/CODE/NK/), which uses some code from Terry Jones 
 (https://github.com/terrycojones/nk-landscapes).



 ##### Our nk_walk.cpp history:
 - v.0.0 2010-08-23: original code by Wim Hordijk.
 - v.1.0 2019-09-12: first versions with the inversions mutations operation.
                   Done by Paul Banse (paul.banse@inria.fr) and
                   Leo Trujillo (leonardo.trujillo@inria.fr)
- v.1.1 2020-02-28: final version to perfom exhaustive exploratory expirical 
                   analysis. Assembled to the "**AWNK** numerical project" (now **SimNK**), with
                   the scripts to run the simulations and data analysis.
                   Done by Leo.
- v.2.0 2020-08-06: calibrations of the numerical experiments and benchmarks. 
                   Exhaustive empirical analysis (first report 2020-10-15).
- v.3.0 2020-11-24: some modifications and the correction of the bias towards 
                   inversions of size two. Done by
                   Guillaume Beslon (guillaume.beslon@inria.fr)
- v.3.1 2021-01-15: benchmarks and test of new parameters. Reorganization
                   of the code and standarization. Incorporation to the **NK_hikes** project. 
- v.3.2 2021-02-04: implements the option to read an initial genome from a 
                   file containing a genome string of size N (for a binary 
                   alphabet A = {0,1}).

