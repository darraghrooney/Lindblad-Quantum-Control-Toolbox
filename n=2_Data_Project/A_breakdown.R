# This function takes a positive-definite 3x3 complex matrix, and decomposes it
# into three parts: the sorted vector of eigenvalues of Re(A), the eigenvector
# matrix of Re(A), and the b-vector representing the imaginary parts of the
# off-diagonal elements.

A_breakdown <- function(A){

  eiginfo = eigen(Re(A), symmetric=TRUE)
  a = eiginfo$values
  rot = eiginfo$vectors*det(eiginfo$vectors)
  b = 2 * t(rot) %*% matrix(c(Im(A)[2,3],Im(A)[3,1],Im(A)[1,2]))

  return(list(a=a,b=b,rot=rot))
}

