Octave/MATLAB code for calculating steady-states and constraints of Lindblad systems
------------------------------------------------------------------------------------

Author: Darragh Patrick Rooney

A quantum system under Lindblad dissipation is contracting on the space of lambda, the vector of density matrix
eigenvalues. lambda obeys a linear ODE: `d(lambda)/dt = Omega*lambda`. The matrix Omega has
off-diagonal elements `w(i,j)` and off-diagonal elements `-sum(i!=j) w(i,j)`. The columns sum 
to zero, so that the sum of lambda remains one. 

Let G be the directed graph where each edge corresponds to a non-zero w(i,j) (pointing from node j to node i). Assume G is connected (if not, each component can be dealt with individually). It is a theorem that if the number of sinks (nodes with zero outgoing edges) is zero or one, the system shrinks to a unique asymptotic state lambda\_inf. Define: `J_k = sum( product( w(i,j) ))` where each product represents a connected tree rooted at node k, and we sum over all such trees. Then the kth component of `lambda_inf` is `J_k / sum_l(J_l)`.

If there are two or more sinks (let m1 = sink number, m2 = source number), there is no unique asymptotic state, but we can say that
1. the components of lambda corresponding to the sources approach zero asymptotically.
2. lambda obeys m2 constraints (the sum(lambda) = 1 constraint is contained therein)

The jth constraint coefficients can be written:

c_j = sum(l) K(j,l)
c_k = 0, 1 <= k <= m1, k !=j
c_l = K(j,l), k > m1

where K(j,l) is defined similarly to J_k, but now the sum is done over tree products. A tree product in this case involves first assigning each source to a sink, then using one rooted tree for each source to write down the product. K(j,l) is the sum of all such tree products such that source l is in a tree rooted at sink j.

### List of .m files:

* `rand_w.m`

This generates a random matrix of rates w(i,j). Each rate is uniform on [0,2] so that the products don't get crazy. A second input can be used to set the first m1 components to be sinks.

* `w_get.m'

Calculates w-rates for a given set of Lindblad operators

* `Omega.m`

This calculates the matrix Omega defined above from a given set of rates.

* `stat_states.m` 

Given a matrix of rates, calculates the asymptotic lambda defined above.

* `rt_w.m` 

This function calculates the sum of products over all rooted trees, where the root and the set of descendent vertices are given as input. 

* `inc_nd.m` 

A short helper function that increments a number on the space (Z_a)^b. Useful as some loops in other functions must be done in such space.

* `multisink.m` 

Calculates the constraint coefficients defined above. From a given m1, m2, and rate matrix, returns the m1 x (m1+m2) matrix of coefficients.

* `constraint_check.m` 

Checks that the constraint matrix C given by `multisink.m` in fact obeys the relation `C^T Omega = 0`.

* `create_3x3.m'

Creates a set of Lindblad operators for a n=6, m1=3 system. The Lindblad operators are jump operators, with the option of having one de-phasing op.

* `upper_3x3.m'

Maps the possible stationary-states when two of the non-sink states are rotated. Stat-states for random unitary (in the non-sink sector only) flags are calculated and plotted for comparison.

* `lower_3x3.m'

Maps the possible stationary-states when two of the sink states are rotated. Stat-states for random unitary (in the sink sector only) flags are calculated and plotted for comparison.

* `flag_3dgen.m'

Generates a random n=3 flag (helper function for the '***3x3.m' files)

* `order_3drates.m'

Permutes indices in a n=3 system so that the rates are ordered (helper function for the '***3x3.m' files).