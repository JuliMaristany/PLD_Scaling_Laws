# Scaling Laws for PLD domains: MD Simulations

## Source code for "Decoding Phase Separation of Prion-Like Domains through Data-Driven Scaling Laws"

We are delighted to share our code with the community. Please use it freely and cite our paper (Preprint: https://doi.org/10.1101/2023.06.14.543914). We are happy to answer any questions and comments by email (jerellejoseph@princeton.edu), and welcome contributions for any updates.

## System requirements

Linux with C++ compilers with MPI. Tested on: CSD3 cclake cluster (https://docs.hpc.cam.ac.uk/hpc/user-guide/cclake.html) with Intel compilers

LAMMPS is required for running MD sims. To build follow instructions on https://docs.lammps.org/Build.html, specifically for the stable version of September 2021. The typical installation time for a desktop computer is roughly 1 hour.

## Demo

To run an example simulation,

1. Move into the demo folder
2. Run lammps using at least 8 cores 

 > mpirun -np 8 ./lmp -in lammps.in

These simulations generate slab trajectories for a particular example PLD (TIA1 Wild Type).

## Reproduction of results: Simulation and analysis of an MD slab trajectory

To reproduce a particular binodal, or to create a new one,

1. Prepare input scripts:
  - Move into the Code/Inputs folder
  - Replace the sequences in the file with the sequence of your desired PLD
  - Run script to obtain lammps config file
2. Run MD simulation:
 - The lammps scripts provided in Code/Simulation take a single PLD as an input and compute the slab geometry from scratch, first by replicating the protein and then by entending the simulation box in one direction. The only user input needed is the name of the config file, and the range of temperatures wished to simulate
3. Analysis scripts are provided in Code/Analysis, to obtain first the density profiles and then the bindals fitted witht he law of rectilinear diameter. More detailed instructions are found in the README.md files on Code/Analysis.

