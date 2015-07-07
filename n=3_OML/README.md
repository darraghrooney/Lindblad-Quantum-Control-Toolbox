Octave/MATLAB code for analyzing n=3 Lindblad Quantum Control systems
---------------------------------------------------------------------

Author: Darragh Patrick Rooney

An n=3 Lindblad quantum control system with fast unbounded control can be characterized by the ODE:

d(lambda)/dt = A'(w)*lambda

where lambda is the 3-vector of eigenvalues, summing to one. w is a matrix of Lindblad rates:

w_(i,j) = sum_k L_(k,ij)

where L_(k,ij) is the ij'th element of the k'th Lindblad operator (which is basis-dependent). A'(w) is a matrix thats depend linearly on the off-diagonal elements of w. The element-sum of lambda is preserved (convention is one), so that the system is actually 2-dimensional. Furthermore, each element must be non-negative, so that our state space is a 2-simplex (technically the quotient of a 3-simplex with the action of S3). This simplex is called T.

This code ignores the quotient technicality and projects the system onto a 2-simplex in R^2. The ODE becomes:

dx/dt = b(w) - A(w)*x

where x is lambda projected, b and A are a vector and matrix respectively, that depend linearly on the off-diagonals of w. The projection of T is called bar{T}. 

The six (off-diagonal) elements of w depend on the basis or flag. If one has complete control over the Hamiltonian, this flag can be viewed as a control variable. An important question is what points on bar{T} can be reached using a set of given flags. Equivalently, if you know the set of flag available, you can calculate in principle the set of w's. If the set is discrete this is practical, but the image of a continuous set (for example, all flags) is not at all straightforward. So a set of discrete flags can be selected for the sake of practicality.

Given a discrete set of control flags, one can plot the arcs that contribute to the boundary of the strong local controllability (SLC) region, or AE3. There are three plotting functions here that calculate the arcs, the resulting boundary and the fill. The globally controllable region is larger, and requires finding critical points on the boundary and iterating the region. Plotting functions are here for only one iteration. More than one was too messy and not of interest. There are additional tools here for plotting trajectories, calculating flag derivatives and visualizing the available tangent vectors. Here follows a complete list. 

### List of .m files:

* `rand_Ls.m`, `rand_GKS.m` and `rand_w.m` are all used for generating randomized system data. The first produces a specified number of Lindblad operators, the second a GKS matrix and the third a w-vector. `flag_gen.m` is a function for producing a random flag.

* `GKStoL.m`, `GKStow.m`, `LtoGKS.m` and `Ltow.m` are conversion functions for converting between GKS, the Lindblad ops, and the w-vector. `principal_w.m` calculates the w-vector for a set of Lindblad operators and the identity flag. `order_rates.m` rearranges a w-vector so that it's asymptotic state is in the primary sector.

* `A_b_rinf.m` calculates A(w), b(w) and rinf = A(w)\b(w) for a given w-vector.

* `AE3_arcs.m` calculates plots the arcs in bar{T} that yield the SLC region. By default, they are plotted for the permutations of the identity flag. Additionally, arcs for randomized flags, along with their permutations, can be plotted. `AE3_fill.m` is similar, but is for plotting the filled regions. It has the option to only show the folded SLC region. `AE3_boundary.m` instead of calculating the arcs, goes ahead and plots the boundary (erasing interior arcs). `arenas.m`and is a graphical helper function that plots the boundary and interior edges of bar{T} around the arcs or fill.

* `BO3_arcs.m`, `B03_boundary.m` and `B03_fill.m` are like the AE3 functions, but are for finding the first generation of higher controllability, i.e. the points that can be attained from AE3 using constant-controls. `boundary_check.m` is a function that checks to see whether this first generation is the last or whether one can escape it.

* `disconnected.m` is a function for plotting the SLC region when A(w) is degenerate. i.e. the graph of non-zero w's is disconnected.

* `DF_w.m` and `DF_Ab.m` calculate the flag derivatives of w, A and b. `DF2_w.m` and `DF2_v.m` calculate the second derivatives of w and v=b-Ax.

* `trajec.m` calculates a trajectory through bar{T} given an initial position and system data. `trajectories.m` will plot multiple trajectories starting from grid points on the boundary.

* `isconemax.m` takes a w-vector and a point in bar{T} and decides whether the resulting six v-vectors form a complete cone.

* `tangents.m` discretizes the flag space, and plots a v-vector for each flag. `tangents_with_princ.m` does the same but adds in the principal v-vectors.

* Helper functions: `flag_inc.m` is an increment function for sweeping through the space of flags. `gPaulis.m` generates the eight generalized Pauli matrices, while `su3off_basis` only generated the six non-diagonal ones. `S3gen.m` generates a rep of the S3 group onto 2x2 matrices, i.e. the linear transformations under which bar{T} is invariant. `rot_curves.m` take a trajectory through bar{T} and calculates the five images of the curve under the S3-action. 
	
### Gallery folder:

The Gallery contains examples of the plotting functions for seven different Lindblad systems. The system data is stored in `Example_data.txt`. For systems 1-5, the result of `AE3_arcs.m`, `AE3_fill.m`, `trajectories.m` are plotted, as well as a combo graph of the AE3 boundary and B03 arcs. For system 1, three examples of `tangents_with_princ.m` are given, one example each is given for system 4 and 5. Systems 6 and 7 are disconnected and show the result of `disconnected.m`. 
