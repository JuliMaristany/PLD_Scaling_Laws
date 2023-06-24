# Scripts to run slab MD simulation using lammps

Input scripts of the format SEQUENCE_i.config should each be present on folders of the same name (SEQUENCE_i/). The lammps scripts present here are parallelized, and can therefore run several slab simulations in parallel, depending on cluster limitations. An example cluster submission script (for CSD3 skylake highmem nodes) is found also in the current folder. 

The lammps command to run the simulations with several sequences in parallel is

> mpirun -ppn $mpi_tasks_per_node -np $np lmp_mpi -partition MxN -i lammps.input

where $mpi_tasks_per_node is the number of tasks per node, $np is the number of processes, M is the number of sequences and N the number of nodes assigned to each sequence simulation.
