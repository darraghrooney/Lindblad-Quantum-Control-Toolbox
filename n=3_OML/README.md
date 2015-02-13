Octave/MATLAB code for analyzing n=3 Lindblad Quantum Control systems
---------------------------------------------------------------------

Author: Darragh Patrick Rooney

An n=3 Lindblad quantum control system with fast unbounded control can be characterized by the ODE:

d(lambda)/dt = A'(w)*lambda

where lambda is the 3-vector of eigenvalues, summing to one. w is a matrix of Lindblad rates:

w_(i,j) = sum_k L_(k,ij)

where L_(k,ij) is the ij'th element of the k'th Lindblad operator (obv. basis-dependent). A'(w) is a matrix thats depend linearly on the off-diagonal elements of w. The element-sum of lambda is preserved (convention is one), so that the system is actually 2-dimensional. Furthermore, each element must be non-negative, so that our state space is a 2-simplex (technically the quotient of a 3-simplex with the action of S3).

This code ignores the quotient technicality and projects the system onto a 2-simplex in R^2. The ODE becomes:

dx/dt = b(w) - A(w)*x

where x is lambda projected, b and A are a vector and matrix respectively, that depend linearly on the off-diagonals of w.


### List of .m files:


* `AE3.m`

	
### Examples folder:

Examples
