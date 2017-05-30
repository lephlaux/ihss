# Inexact hierarchical scale separation

*Hierarchical scale separation* (HSS) is an iterative two-scale approximation method for large sparse systems of linear equations arising from *discontinuous Galerkin* (DG) discretizations. HSS splits the linear system into a coarse-scale system of reduced size corresponding to the local mean values of the solution, and a set of *decoupled* local fine-scale systems corresponding to the higher order solution components. This scheme then alternates between coarse-scale and fine-scale system solves until both components converge. The motivation of HSS is to promote parallelism by decoupling the fine-scale systems, and to reduce the communication overhead from classical linear solvers by only applying them to the coarse-scale system.

We propose a modified HSS scheme (*"inexact HSS", "IHSS"*) that exploits the highly parallel fine-scale solver more extensively and only approximates the coarse-scale solution in every iteration thus resulting in a significant speedup. The tolerance of the coarse-scale solver is adapted in every IHSS cycle, controlled by the residual norm of the fine-scale system. Anderson acceleration is employed in the repeated solving of the fine-scale system to stabilize the scheme.

## Third party code
The solver function `ihss.m` requires an implementation of the Anderson acceleration. We recommend Homer Walker's `AndAcc.m`, which can be found [here](https://users.wpi.edu/~walker/Papers/anderson_accn_algs_imps.pdf). 

## Documentation
Run `publish('ihss.m')` in the MATLAB terminal.

## Contact
florian dot frank at rice dot edu
