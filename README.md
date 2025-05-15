# Master's Thesis Project
## Preconditioning strategies for simulating very large floating structures

In this repository, I have added specific code that I have used to generate results in the Master's Thesis report.

There are four sections - 
1. Overlapping methods.
2. Non-overlapping methods.
3. Fluid Structure Interaction Toy problem serial results.
4. Fluid Structure Interaction Toy problem parallel results.
5. Simulation of the floating membrane.

To use this code, first install the Trilinos library as a shared library. You then need to update the paths in the cmakelists
file and some bash scripts. You need to ensure that Julia (MPI.jl) and Trilinos are using the same MPI library.