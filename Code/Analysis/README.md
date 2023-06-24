# Analysis scripts to obtain critical temperatures from a slab MD simulation

- tail.sh is a bash script that grabs the final 500 ns of the simulation, and splits it in chunks of 100 ns. For the simulation lenghts we have,  we have proved that this yields a well-converged result, but for bigger systems, these numbers may need to be altered
- ppd.sh calls on calculate.
