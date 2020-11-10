# lattice-percolation
Lattice percolation in C using MPI and OpenMP for parallel computing.

## How to run
1. Download percolation.c
2. Enter `mpicc -o percolation percolation.c -fopenmp` to compile percolation.c and create and executable file percolation.
3. Enter `mpirun -np 3 ./percolation` to execute the application.

