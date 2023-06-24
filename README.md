# Universal Codes for PLD domains: MD Simulations

## Source code for "Universal Predictive Scaling Laws of Phase Separation of Prion-Like Low Complexity Domains"

We are delighted to share our code with the community. Please use it freely and cite our paper (Preprint: https://doi.org/10.1101/2023.06.14.543914). We are happy to answer any questions and comments by email (jerellejoseph@princeton.edu), and welcome contributions for any updates.

## System requirements

Linux with C++ compilers with MPI. Tested on: CSD3 cclake cluster (https://docs.hpc.cam.ac.uk/hpc/user-guide/cclake.html) with Intel compilers

LAMMPS is required for running MD sims. To build follow instructions on https://docs.lammps.org/Build.html, specifically for the stable version of September 2021.

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
  - 
2.  

