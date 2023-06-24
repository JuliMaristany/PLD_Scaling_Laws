# Analysis scripts to obtain critical temperatures from a slab MD simulation

- tail.sh is a bash script that grabs the final 500 ns of the simulation, and splits it in chunks of 100 ns. For the simulation lenghts we have,  we have proved that this yields a well-converged result, but for bigger systems, these numbers may need to be altered
- calculate.py computes the mean density profile of all chunks
- ppd.sh calls on denfitting.py, which fits the density profile to obtain the higher and lower density values (for the condensed and dilute phase respectively). The names of the files need to be updated to your naming convention.
- binodal.py plots the binodals and computes critical temperature
