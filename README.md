# Inexact hierarchical scale separation

*Hierarchical scale separation* (HSS) is an iterative two-scale approximation method for large sparse systems of linear equations arising from *discontinuous Galerkin* (DG) discretizations. HSS splits the linear system into a coarse-scale system of reduced size corresponding to the local mean values of the solution, and a set of *decoupled* local fine-scale systems corresponding to the higher order solution components. This scheme then alternates between coarse-scale and fine-scale system solves until both components converge. The motivation of HSS is to promote parallelism by decoupling the fine-scale systems, and to reduce the communication overhead from classical linear solvers by only applying them to the coarse-scale system.

We propose a modified HSS scheme [1] (*"inexact HSS", "IHSS"*) that exploits the highly parallel fine-scale solver more extensively and only approximates the coarse-scale solution in every iteration thus resulting in a significant speedup. The tolerance of the coarse-scale solver is adapted in every IHSS cycle, controlled by the residual norm of the fine-scale system. Anderson acceleration is employed in the repeated solving of the fine-scale system to stabilize the scheme.

## Third party code
IHSS requires Anderson acceleration for the repeated solving of the fine-scale problems within one HSS cycle.  The implementation `AndAcc.m` [2,3] was uploaded into our repository by kind permission of the author Homer F. Walker.

## Documentation
Run `publish('ihss.m')` in the MATLAB terminal.

## Contact
florian dot frank at rice dot edu

## References
[1] [C Thiele, M Araya-Polo, FO Alpak, B Rivière, F Frank, Inexact hierarchical scale separation: A two-scale approach for linear systems from discontinuous Galerkin discretizations, Computers and Mathematics with Applications, 2017.](http://dx.doi.org/10.1016/j.camwa.2017.06.025)

[2] [HF Walker, Anderson Acceleration: Algorithms and Implementations, 2011](https://users.wpi.edu/~walker/Papers/anderson_accn_algs_imps.pdf)

[3] [HF Walker, P Ni, Anderson Acceleration for Fixed-Point Iterations, SIAM J. Numer. Anal., 49(4), 1715–1735, 2011.](https://doi.org/10.1137/10078356X)
