Octave/MATLAB code for analyzing n=4 Lindblad Quantum Control systems
---------------------------------------------------------------------

Author: Darragh Patrick Rooney

An n=4 Lindblad quantum control system with fast unbounded control can be characterized by the ODE:

d(lambda)/dt = Omega(w)*lambda

where lambda is the 4-vector of eigenvalues, summing to one. w is a matrix of Lindblad rates:

w_(i,j) = sum_k L_(k,ij)

where L_(k,ij) is the ij'th element of the k'th Lindblad operator (obv. basis-dependent). Omega(w) is a matrix that depends linearly on the off-diagonal elements of w. The element-sum of lambda is preserved (convention is one), so that the system is actually 3-dimensional. Furthermore, each element must be non-negative, so that our state space is a 3-simplex (technically the quotient of a 3-simplex with the action of S4).

This code ignores the quotient technicality and projects the system onto a 3-simplex in R^3. The ODE becomes:

dx/dt = b(w) - A(w)*x

where x is lambda projected, b and A are a vector and matrix respectively, that depend linearly on the off-diagonals of w.

The goals of this code are essentially (1) to plot the STLC set of our system, if one views the control set as the 24 permutations of the basis set, and (2) for standard Lindblad ops,  to compare stat. orbits of random flags with those of the one-parameter rotations. 

### List of .m files:

* `AE4.m`

	This is the main function, responsible for plotting the *likely* boundary of the STLC set. The actual boundary is composed of 2024 "curved triangles" (i.e. quadratic images of 2-simplices). This number of meshes can overwhelm the CPU - an alternate strategy is to choose a subset of 92 surfaces corresponding to neighboring triples (instead of all triples on Z_24). Usually, but not always, the remaining surfaces remain interior to this boundary, so this function plots only these 92 surfaces. Usually a choice of grid size = 20 is sufficient to get a nice-looking, manageable plot on a cheap modern laptop.

* `AE_alt.m`

	This "alt" version (attempts to) plot all 2024 surfaces. No example plots included here were generated with this function. Included only for completeness. Use at own risk.

* `AE_bound.m`

	To check that the `AE4.m` plot actually characterizes the STLC set, one can plot only the edges of the 2024 surfaces (there are 276 such edges). This function does that. It can be combined with the plot in `AE4.m` to overlay the abridged plot with the unabridged edge-plot. The Examples folder contains one example where this is useful.

* `rand_w.m`

	Produces a randomized w-matrix. Diagonals are zero. Off-diagonals are uniform over the interval [0,2]. 

* `rand_stdLs.m`

	Produces randomized standard Lindblad ops. i.e. 12 different jump operators, and one de-phasing operator.

* `Ltow.m`

	Converts Lindblad operators to a w-matrix given a flag.

* `flag_gen.m`

	Generates a flag on C^4. Either takes 12 input parameters, or selects them randomly.

* `boundaries.m`

	This function plots the edges of the 4-simplex corresponding to the boundary of our space. The vertices are:
	* (1/sqrt(2), 1/sqrt(6), 1/sqrt(12))
	* (-1/sqrt(2), 1/sqrt(6), 1/sqrt(12))
	* (0, -2/sqrt(6), 1/sqrt(12))
	* (0, 0, -3/sqrt(12))

* `A_b_rinf.m`

	This function calculates A(w) and b(w). Also returns the asymptotic state `rinf = A\b`.

* `order_rates.m`

	It is often useful to permute the Lindblad rates so that the asymptotic state falls into the primary sector (corresponding to descending eigenvalues). This function does that and returns the input rates J_k.

* `rot_surf.m`

	Given a surface in which one of the vertices is 1, this function transforms this surface according to one of the S4-generating matrices (which one is specified as input). This function is required for `AE4.m`.

* `S4gen.m`

	This function returns the 24 matrices for the representation of S4 on R_3. These 24 matrices permute the 4 vertices of the simplex drawn in `boundaries.m`. 

* `std_bound_test.m`

	This function compares the stationary orbits of the rotational flags to any number of randomized flags. This is intended for jump and de-phasing Lindblad ops only. Only the primary sector is plotted. In that sector, there are six curves corresponding to six pairs of basis vectors. The number of random flags is taken as an input parameter.

### Examples folder:

Examples of plots for `AE4.m`, `AE4_bound.m` and `std_bound_test.m` are included. There are eight examples of the first two and four of the latter. For the second example, the "skeleton" contains pieces that are external to the `AE4.m` surfaces, and so the skeleton and surfaces are super-imposed on a third plot for that example. This does not apply to the other seven plots of `AE4.m`. The specific w-matrices/L-parameters corresponding to all plots are stored in the text file `example_data.txt`.
