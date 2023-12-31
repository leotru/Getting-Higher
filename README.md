# Getting Higher

This repository contains supporting data and codes for the paper

### Getting higher on rugged landscapes: Inversion mutations open access to fitter adaptive peaks in NK fitness landscapes
by L. Trujillo, P. Banse and G. Beslon 

PLoS Computational Biology, October 31, 2022

https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010647

#### Abstract
Molecular evolution is often conceptualised as adaptive walks on rugged fitness
landscapes, driven by mutations and constrained by incremental fitness selection. It is
well known that epistasis shapes the ruggedness of the landscape’s surface outlining
their topography (with high-fitness peaks separated by valleys of low fitness genotypes).
However, within the strong selection weak mutation (SSWM) limit, once an adaptive
walk reaches a local peak, natural selection restricts passage through downstream paths
and hampers any possibility of reaching higher fitness values. Here, in addition to the
widely used point mutations, we introduce a minimal model of sequence inversions to
simulate adaptive walks. We use the well known NK model to instantiate rugged
landscapes. We show that adaptive walks can reach higher fitness values through
inversion mutations, which, compared to point mutations, implies that the evolutionary
process escapes local fitness peaks. To elucidate the effects of this chromosomal
rearrangement, we use a graph-theoretic representation of accessible mutants and show
how new evolutionary paths are uncovered. The present model suggests a simple
mechanistic rationale to analyse escapes from local fitness peaks in molecular evolution
driven by (intragenic) structural inversions and reveals some consequences of the limits
of point mutations for simulations of molecular evolution.

The preprint is available on (soon on arXiv)

##### Repository content
Codes:

nk_walk.cpp: Main routine to simulate adative walks on the Stuart Kauffman's NK-fitness landscape model, with point mutations (as usual) and inversions (as new).

mean_fitness_full_plot.py: Python code to calculate and plot the mean fitness.

full_mutational_network.m: Matlab scritp to calculate mutational networks.

fitnessGraphs.m: Matlab script to calculate the fitness networks.

Data: The data used to plot figures 3 and 4 of the manuscript. 

Figures: Figures 3 and 4 of the manuscript (also two supplementary figures of the mean fitness for N=50 and 80). 






