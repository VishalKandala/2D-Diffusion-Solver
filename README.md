# PETSCHeat
In this work, the 2-D Diffusion equation is solved on a structured grid using the Finite
Difference Method, This is done using a parallelized code, utilizing the PETSc library’s C
interface,run on the Terra cluster. Krylov Subspace Solvers were utilized and the solutions
obtained (Time to convergence) through Gauss Siedel Iteration and the standard GMRES

method are compared, along with a strong scaling analysis of both those solutions is illus-
trated.

1 Introduction
The objective of this work is to solve the 2-D Diffusion equation on a structured grid using the
Parallel PETSc Library, this section outlines the problem domain, boundary conditions and the
anatomy of the governing differential equation.
1.1 Domain and Boundary Conditions
A rectangular domain of unit thickness, with length 0.4 m and width 0.3m is bounded by dirichlet
boundary conditions on all four sides, the domain can be seen in Fig.1. The material is assumed
to have thermal diffusivity (α) of 11.234 x 10−5m2/s and a thermal conductivity of 280 W/m.K,
the plate is assumed to be homogenous and as having constant properties across the domain.

Figure 1: The Computational Domain

1

The Dirichlet boundary conditions imposed on the system mean that the boundaries are all
held at temperatures shown in Table.1, where the edge which is held at the label can be seen in
Fig.1.

Label Temperature
T1 313 K
T2 273 K
T3 283 K
T4 273 K
Table 1: Boundary Temperatures

1.2 Governing Differential Equation
The 2-D Steady State Heat Equation or the diffusion equation is shown in Eq.1, where T is the
dependent variable and x,y are the independent variables.

∂
2T
∂x2
+
∂
2T
∂y2
= 0 (1)
This is a linear,homogenous,Second Order Elliptic Partial Differential Equation (B2 − 4AC < 0)
and hence exhibits two-way street behaviour or Boundary Value Problem (BVP) behaviour. This
means that the information propagation happens at infinite speed and in all directions.
Diffusion phenomenon involves a gradual reduction of gradients across the domain and hence

results in an averaging or smoothening effect and that is what is expected from the solution irre-
spective of the methodology implemented.

1.3 PETSc Linear Algebra Solvers
PETSc[1] is a Scientific Computing library that contains highly optimized, parallelized routines for
Matrix and Vector Operations through custom objects (Vec,Mat,KSP) as well as the Distributed
Array (DM) data structure which links the linear algebra solvers to the matrix and vector objects.
2 Methodology
Numerical solution of a Partial Differential Equation with non-homogenous boundary conditions
can be obtained through the following three steps, each of which will be detailed in this section
1. Grid Generation.
2. Obtaining a System of Linear Equations.
3. Solving a system of Linear Equations.
2.1 Grid Generation
A structured, uniform, node based grid with a star stencil has been implemented to solve this
problem, in the program, the variable i has been assigned to grid points along x-axis and the
variable j has been assigned to grid points along y-axis.The stencil can be seen in Fig.3 and a 8x6
grid can be seen in Fig2. The grid spacing ∆x and ∆y for an N x M grid are given in Eq.2

∆x =
L
N − 1
; ∆y =
W
M − 1

(2)

2

Figure 2: The Finite Difference grid utilized to solve the diffusion equation with the boundary
grid points marked in green and the internal grid points marked in red, along with their indices.

Figure 3: The Finite Difference star stencil implemented in this project

3

2.2 Finite Difference Scheme
To generate the linear system of equations, a 2-D second order central differencing scheme was
implemented along both the x and y axes, on the star stencil, as shown in Eq.4

Letz = f(x, y) (3)
∂
2
z
∂x2
=
z(i + 1, j) + z(i − 1, j) − 2z(i, j)

∆x
2

(4)
To derive the accurate finite difference equations at various locations, the grid has been divided
into 10 different regions which are illustrated in Fig.4, the finite Difference equations obtained in
each region are described in the following sub sections.

Figure 4: Boundary Condition Implementation through splitting the domain into multiple regions
2.2.1 Internal Nodes
For internal nodes i.e i>1,j>1 and i<N-2,j<N-2, discretizing the diffusion equation using the
second order central differencing approach shown above leads to linear equations shown in Eq.7

1
∆x
2
= Dx (5)
1
∆y
2
= Dy (6)
2(Dx + Dy)T(i, j) − DxT(i − 1, j) − DxT(i + 1, j) − DyT(i, j + 1) − DyT(i, j − 1) = 0 (7)
2.2.2 Boundary Nodes
Since Dirichlet Boundary conditions are to be applied, the finite difference equation obtained at
the boundaries i.e i=0 or j = 0 or i=M-1 or j=N-1 simply places the boundary temperature in it’s
location as shown in Eq.8, Eq.9,Eq.10 and Eq.11.

T(i, 0) = T1(i > 0;i < M − 1) (8)
T(0, j) = T2(j > 0; j < N − 1) (9)
T(i, N − 1) = T3(i > 0;i < M − 1) (10)
T(M − 1, j) = T4(j > 0; j < N − 1) (11)

4

At the corners i.e the singularities, these points are not solved for anyway and hence, to maintain
numerical consistency, decoration values obtained by averaging the boundary temperatures of the
edges that coincide in that corner are placed.
2.2.3 Regions 1,2,3,4
These regions correspond to the nodes that are close to the edges at the left,right,top and bottom

boundaries respectively, here, one of the temperatures is a known value and is equal to the bound-
ary value, hence, these equations have a constant term and three unknowns with co-efficients. The

equations for Regions 1,2,3 and 4 are shown in Eq.12,Eq.13,Eq.14 and Eq.15 respectively.
2(Dx + Dy)T(1, j) − DxT(2, j) − DyT(1, j + 1) − DyT(1, j − 1) = DxT2
(12)
2(Dx + Dy)T(M − 2, j) − DxT(M − 3, j) − DyT(M − 2, j + 1) − DyT(M − 2, j − 1) = DxT4
(13)
2(Dx + Dy)T(i, N − 2) − DyT(i, N − 3) − DxT(i + 1, N − 2) − DxT(i − 1, N − 2) = DyT3
(14)
2(Dx + Dy)T(i, 1) − DyT(i, 2) − DxT(i + 1, 1) − DxT(i − 1, 1) = DyT1
(15)

2.2.4 Regions 5,6,7,8
These nodes are near corners of the domain and hence, two of the neighbouring nodes would have
known (boundary) temperature values, which leads to Eq.16,Eq.17,Eq.18 and Eq.19 respectively.
2(Dx + Dy)T(1, 1) − DxT(2, 1) − DyT(1, 2) = DyT1 + DxT2
(16)
2(Dx + Dy)T(M − 2, N − 2) − DxT(M − 3, N − 2) − DyT(M − 2, N − 3) = DyT3 + DxT4
(17)
2(Dx + Dy)T(M − 2, 1) − DxT(M − 3, 1) − DyT(M − 2, 2) = DyT1 + DxT4
(18)
2(Dx + Dy)T(1, N − 2) − DxT(2, N − 2) − DyT(1, N − 3) = DyT3 + DxT2
(19)

2.3 Solving System of Linear Equations
Formulating as shown above would lead to a system of linear equations of the form Ax=b which

can be solved either directly (Gauss Elimination) or through iterative methods (Jacobi,Gauss-
Seidel etc).

More often than not, unless the system is a very simple system with very low condition number
i.e the ratio of the highest and lowest singular values of the co-efficient matrix, for example, a
Tri-Diagonal system, these systems of linear equation are not amenable to direct solution and
hence iterative solvers are preffered.

Some methods such as the Strongly Implicit Procedure combines both direct and indirect (itera-
tive) methods to obtain solutions. In this section, the Gauss-Siedel iterative method is discussed

and it’s implementation using the PETSC linear solver class ”KSP”[2] is detailed, there will also
be a brief summary of Krylov Sub-Space methods and the GMRES[6] method that is implemented
as the default solver in PETSc.
2.3.1 Gauss Siedel Iteration
This involves iteratively minimizing the norm of the residual for each grid point, consider an
internal point i,j, then the linear equation at this point is shown in Eq.20,here we observe that
the solution ”propogates” to the east and north and subsequently as the solution already was
undertaken in the current iteration at the south and west neighbour of the current point, updated

5

values are available there and the values from previous iteration are available at the north and
east neighbours. This is also illustrated in Fig.5
T
k+1(i, j) = DxT

k+1(i − 1, j) + DxT
k
(i + 1, j) + DyT

k+1(i, j − 1) + DyT
k
(i, j + 1)

2(Dx + Dy)

(20)

Figure 5: Information propogation in a Star-Stencil centered around the point i,j when Gauss-
Siedel Iteration is used to solve the system

When considering the whole system, this can be represented as the splitting or ”precondition-
ing” of the co-efficient matrix.

Consider the simplest solver i.e the richardson solver which is represented in Eq.22, which simply
adds the residual to the current iteration to obtain the next iteration and minimizes the residual,
this has been modified by adding a Pre-conditioner P as shown in Eq.23,accordingly, there are
many ways in which a linear system can be factorized/split or transformed to better condition for
solvability and an array of options are available in PETSC through the PC objects.
Parallelizing Gauss-Seidel iteration could be quite challenging and the most direct methods are
the wavefront method and coloring method[4], however, it can be seen as the following split of
the system shown in Eq.24, which means, by applying this as a preconditioner to the modified
richardson iteration solution(Eq.23), we can obtain a gauss-siedel solution equivalent[3] to Eq.20
as demonstrated in Eq.25

Ax = b (21)

x
k+1 = x
k + (b − Ax
k
) (22)

x
k+1 = x
k + P(b − Ax
k
) (23)
A = L∗ + U (24)
x
k+1 = L
−1
∗
(b − Ux
k
) (25)

6

A slightly better iteration method which converges quicker is known as Successive Over/Under

Relaxation and it simply multiplies a constant ω to the pre-conditioned residual, if this ω is cho-
sen to be 1, then the solver is identically Gauss-Siedel and that is the way in which Gauss Siedel

Iteration can be implemented in Petsc by using ”richardson” solver, ”sor” pre-conditioner with
”omega” set to one and non-zero initialization that accelerates convergence.
Another way to implement Gauss-Seidel in PETSc is to use the ”preonly” option available with
the KSP solver object and using an ”sor” preconditioner and setting ”omega” to one as before.
2.3.2 Krylov Sub-Space Methods

Basic Iterative solution methods such as Jacobi, Gauss-Seidel and SOR are only valid for diag-
onally dominant systems and the rates of convergence can be quite low, hence more advanced

methods are necessary to solve systems which are not as ”nice”.
Consider the initital solution x0, the co-efficient matrix A, then the vectors AX0,A2x0,A3x0 upto
Akx0 span the ”krylov Sub Space” of the system and various methods that involve these subspaces
have been proven to be more efficient and adept at solving systems of linear equations, even if
they are singular, assymmetric and/or not diagonally dominant.
The default solver that is implemented in PETSc and one of the most popular Krylov Subspace
solver, one that does not assume any prior knowledge of the system matrix A is the Generalized
Minimum Residual Method (GMRES), this involves minimizing the residual over the entire Krylov
Subspace as shown in Eq.26 and the linear combination that is part of the Krylov Subspace that
minimizes the residual is the solution vector xk as shown in Eq.27

minckA(c0x0 + c1Ax0 + c2A

2x0 + . . . ckA

kx0) − bk (26)

xk = c0x0 + c1Ax0 + c2A

2x0 + . . . ckA

kx0 (27)
The GMRES[6] Solver with incomplete LU Factorization is the default setting for the KSP Solver
object available in the PETSc library,for the scaling analysis, this solver will be compared with
the Gauss-Siedel Method.
3 Results and Discussion
Figure.6 shows the temperature contour of the steady state that the system described by the
given domain would settle in, given that only diffusion is at play. This particular contour plot
is produced with a grid of 400 nodes along the x-axis and 300 nodes along the y-axis that was
solved on 4 Intel Xeon E5 Processor cores on the Terra Computing cluster. The solver used is the
Generalized Minimal Residual Method with Incomplete LU Factorization pre-conditioning
From the figure, it can be observed that the boundaries are maintained at 313K,273K,283K and

273K respectively and the diffusion phenomenon spread the temperature distribution around ef-
fectively, the highest temperature in the domain is at the bottom edge and the least temperature

is at the left and right edges, all the interior temperatures were between these highs and lows
which was as expected. The temperature at the center of the domain is close to 290K and the
temperature contours form a sort of progression from the bottom edge to the top which is what
you would expect the diffusion operator to do.

7

Figure 6: Steady State Temperature Contour obtained with a grid size of 400x300, solved using
GMRES and ILU Preconditioning.
Figure.7 shows the temperature contour obtained with a coarser 100x100 grid implemented
on 4 Intel Xeon E5 Processor Cores using the Gauss-Siedel Iteration Method, and even in this
contour plot the same features can be observed as in Fig.6 i.e the averaging effect of the diffusion
process is apparent.

Figure 7: Steady State Temperature Contour obtained with a grid size of 100x100, solved using
Gauss-Siedel Iteration.

Case N M Processors Solver Preconditioner Tolerance Iterations to Convergence
1 400 300 4 GMRES Incomplete LU 1e-5 456
2 100 100 4 Richardson SOR (ω = 1) 1e-5 3510
Table 2: Detailed Description of the cases depicted in Fig.6 and Fig.7

8

3.1 Scalability Analysis
The algorithm has been parallelized to take advantage of the capabilities of distributed computing
power available through the Terra cluster which houses 28 Intel Xeon E5 cores in each node.
To truly be able to gauge the benefit of parallelization, strong and weak scaling are defined, if
an algorithm demonstrate strong scaling which is a stricter condition that would be ideal, to
understand this metric, we look at a metric called speed up (Eq.28) which is regularly employed
to analyze HPC performance. Here, t1 is the time taken to execute the program by one particular
core/processor where as tN is the time taken to execute the program by N cores/processors.

speed − up =
t1
tN

(28)

Strong scaling is gauged by analysing the variation in speed-up as the number of processors in-
creases, ideally this would be a linear relationship, however, as the number of processors increases

the time taken to communicate between processors increases even as the time each processor takes
to compute it’s share of the code comes down and hence speed-up flattens out,this is known as
Amdahl’s law[5], practically though, performance tapers off even before and even fluctuates due

to constraints placed by the architecture of the HPC system in use and the interconnect (commu-
nication) system in use.

In this particular project, the strong scaling performance of the gauss-seidel algorithm is tested
with the problem size being set at 100x100 grid, and the GMRES Solver is analyzed with a grid
size of 400x300. The computation times and the speed-up trends can be seen in Table.3, Table.4,
Fig.8 and Fig.9 respectively.

N M Processors Execution Time (s)
100 100 1 0.63
100 100 5 0.18
100 100 10 0.13
100 100 15 0.1
100 100 20 0.09
100 100 25 0.06

Table 3: The Problem size,Number of Processors used and the execution times in each case when
Gauss-Siedel Algorithm was implemented to analyze scaling.
It has to be noted that beyond the 100x100 grid size, Gauss-Siedel solver is not converging
even beyond 10,000 i.e the convergence rate is substantially low even though the trend towards
convergence can be observed.

N M Processors Execution Time (s)
400 300 1 1.61
400 300 5 0.37
400 300 10 0.24
400 300 15 0.16
400 300 20 0.14
400 300 25 0.26

Table 4: The Problem size,Number of Processors used and the execution times in each case when
GMRES Algorithm was implemented to analyze scaling.

9

Figure 8: Speed-Up vs No.of Processors when solving the diffusion equation on a 100x100 grid
using the Gauss-Siedel Method

Figure 9: Speed-Up vs No.of Processors when solving the diffusion equation on a 400x300 grid
using GMRES Method

10

It can be observed from both these plots that the highest speed-up acheived by GMRES solver
is much higher than that achieved by the Gauss-Siedel solver and this speed-up is also achieved
with far fewer CPUs, suggesting that the efficiency of GMRES is far better than Gauss-Siedel.

The drop in speed-up observed in Fig.9 can be explained, as mentioned before, by the archi-
tecture of the nodes in the Terra cluster, Each node has two sockets with 14 cores each that share

their L3 Cache, the memory that is closest to the cores and can be accessed very quickly, hence,
upto, around 14 processors, since the transmission overheads are less since all the cores share the
same memory,we observe higher speed-up, however, beyond this, memory needs to be shared and
communicated between the two sockets which together make up one node. In case of gauss-siedel

iteration, the solution time that each core spends has drastically come down, balancing transmis-
sion overheads, resulting in increased speed-up.

Figure.10 and Fig.11 show the variation of the number of iterations it took for the solution to con-
verge with the number of processors/cores involved in the computation, what is observed is that

for Gauss-Siedel iteration, irrespective of the number of processors used, the number of iterations
settles at a value of around 3550 and flattens out, where as in case of GMRES, there is significant
variation that is also consistent with the speed-up trend observed in Fig.9, The no.of iterations
required to converge declines when the maximum compute that shared memory was employed and
this needs to be explored further by the author.

Figure 10: No.Iterations to Convergence Plotted for the Gauss-Siedel solver on a 100x100 grid

11

Figure 11: No.Iterations to Convergence Plotted for the GMRES solver on a 400x300 grid
4 Summary and Conclusions
The 2-D Diffusion equation was solved successfully over large uniform structured grids using the
PETSc library interface with c language, The nature of the 2D Diffusion equation was explored,
a 3 step process involving Grid Generation, Finite Difference Equation generation and solution of
the thus produced linear system of equations was undertaken.

The grid generated was a structured, uniform grid that spanned the entire domain, it was imple-
mented using the Distributed Array (DMDA) object available as part of the PETSc library, a star

stencil ( Five-point stencil) was implemented.
To generate the difference equations, 2nd order Central Differencing scheme was employed, to
implement the boundary conditions, the domain was split into 8 different regions of which regions
1,2,3,4 correspond to the west,east,north and south edges respectively, regions 5,6,7,8 correspond
to the soth-west,north-east,south-east and north-west corners of the domain respectively.
The system of linear equations of the form Ax=b was obtained through the differencing scheme
and the boundary conditions, to solve this, iterative methods were chosen and specifically two
methods i.e Gauss-Siedel iteration which belonged to the ”primitive” iterative schemes and the

Generalized Minimal Residual (GMRES) Method which is a Krylov Subspace method and as-
sumes no prior knowledge of the co-efficient matrix A i.e sparsity, diagonal dominance, postiive

definiteness.
To implement both these solvers, PETSc’s KSP solver object was utilized, to realize the gauss
siedel solver, ksp type had to be set to the richardson solver with SOR preconditioning, which
when the relaxation factor ω is set to unity, is identical to Gauss-Siedel iteration.
The results did show the averaging effect of the diffusion phenomenon and as were expected,
a strong scaling analysis was done for both Gauss-Siedel ( on a 100x100 grid) and GMRES (on
a 400x300 grid) and it was found that GMRES outperforms Gauss-siedel in terms of maximum
speed up as well as efficiency and Gauss-Siedel requires far more ( at least an order of magnitude

12

higher) iterations to converge for the same problem size.
References

[1] Satish Balay et al. “Efficient Management of Parallelism in Object Oriented Numerical Soft-
ware Libraries”. In: Modern Software Tools in Scientific Computing. Ed. by E. Arge, A. M.

Bruaset, and H. P. Langtangen. Birkh ̈auser Press, 1997, pp. 163–202.
[2] Satish Balay et al. PETSc Users Manual. Tech. rep. ANL-95/11 - Revision 3.15. Argonne
National Laboratory, 2021. url: https://www.mcs.anl.gov/petsc.
[3] Satish Balay et al. PETSc Web page. https : / / www . mcs . anl . gov / petsc. 2021. url:
https://www.mcs.anl.gov/petsc.
[4] Hadrien Courtecuisse and J ́er ́emie Allard. “Parallel dense gauss-seidel algorithm on many-core
processors”. In: 2009 11th IEEE International Conference on High Performance Computing
and Communications. IEEE. 2009, pp. 139–147.
[5] John L. Gustafson. “Amdahl’s Law”. In: Encyclopedia of Parallel Computing. Ed. by David
Padua. Boston, MA: Springer US, 2011, pp. 53–60. isbn: 978-0-387-09766-4. doi: 10.1007/
978-0-387-09766-4_77. url: https://doi.org/10.1007/978-0-387-09766-4_77.
[6] Youcef Saad and Martin H Schultz. “GMRES: A generalized minimal residual algorithm for

solving nonsymmetric linear systems”. In: SIAM Journal on scientific and statistical comput-
ing 7.3 (1986), pp. 856–869.

13

Appendices
A Code and Software Development
The entire project was written in C,Python (Post Processing) as well as bash scripts and SLURM
Job files on the Terra cluster.
PETSc library along with MPI were used on the C programs, numpy, system,os and matplotlib
libraries were utilized to develop the post processing code in Python3.
A.1 File System
The main directory is /scratch/vishalkandala/PETSC/2DIFF/ which contains the following
sub-directories and files:
• src which in-turn contains the libs directory inside which are the main source files 2diff.c,
the makefile and the python post-processing file post
• build which contains the executable 2diff.exe which the make file is src/libs/ would remove
and replace whenever make is called.
• data which contains the output files i.e the solution output from the solvers as wel as the
timegmres.csv and timegs.csv which are updated every time a new run is conducted
and the execution time, iterations to convergence,grid sizes (along x and y) and the no.of
processors used are all appended.
• logs is a directory that contains output files from solves that are of three different kinds
namely, dm.txt files which contain DM information for that run, profile.txt that had all
the profiling information related to the run and res.txt which had the residual output after
every iteration.
• batchproc contains two directories jobs which houses the various SBATCH job files used

while the outputs directory houses the output files generated when a job is done/can-
celled/exits with an error.

• plts is a folder that contains all the plots including the contour plots as well as the scaling
and iteration plots.
• sweep is an executable that deletes all data files, all plots and all the output files.

• cleanmake is an executable that changes directory to src/libs/ and calls make, while re-
placing the existing executable in the build directory

• singlerun is an executable that takes in 6 arguments, namely, number of processors, x-grid
size, y-grid size,contour plot choice (0 means do not plot contour, 1 means plot and produce
an eps file and 2 means plot and produce a png file.) and the speed-up plot choice.
• multipost is an executable that loads the matplotlib module and passes the x,y grid values
as well as the choices for plot and speed-up plot to the python program post.
• batchrun is an executable that changes directory to /batchproc/outputs and submits a job
request using the lsf file in /batch/jobs/
• strongscale is an executable that runs the cases discussed in this project on various number
of processors and stores/creates data in the /data/ directory.
