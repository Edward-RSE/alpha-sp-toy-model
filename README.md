# Toy models for Python

This project includes a selection of toy models for Python.

## `num-int`

This toy model has been made to experiment with the numerical integrators in
GSL to compute the spontaneous recombination coefficient for macro atoms in
the MCRT program Python. 

The hope is that we can find a good compromise between accuracy and speed, as 
the current method in use is becoming too expensive.

## `node-share`

This toy model is used to experiment with using remote memory access (RMA)/node 
sharing/shared memory in MPI to reduce the memory footprint of Python in
MPI. The core problem is there are large structures in Python which do not need
to be duplicated for each rank. With RMA, only one rank needs to read in, for 
example, atomic data and all the other ranks can access the data from the 
root rank.
