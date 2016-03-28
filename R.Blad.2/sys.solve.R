# This file contains the followings functions:
#
# calc.x()
# calc.y()
# calc.C()
# calc.C.helper()
# calc.C.helper.helper()

# For given hamiltonians and measured Bloch vectors, this function calculates x,
# the vectorized system parameters

calc.x <- function(h, n){

  # Calculate C and y
  C = calc.C(n)
  y = calc.y(h, n)
  
  # Calculate x from C and y. Allow for the possibly that C is singular
  x = try(solve(t(C)%*%C,t(C)%*%y))
  
  # If C is singular, return zeros for all parameters
  if (class(x) == "try-error"){
    return(rep(0,9))
  }
  # Specify which parameters are which
  row.names(x) = c("b1","b2","b3","a11","a22","a33","a23","a31","a12")
  
  return(x)
}

# This function calculates y, which is basically a stack of cross-products

calc.y <- function(hs, ns){

  # Matrixify hamiltonians and Bloch vectors (else there are errors)
  hs = as.matrix(hs)
  ns = as.matrix(ns)
  y = matrix(0, dim(ns)[1]*3,1)
  for (j in 1:dim(ns)[1]){
    y[(j-1)*3+1:3,] = -cross.mat(hs[j,]) %*% ns[j,]
  }
  return(y)
}

# This function calculate the matrix C, which is 3M times 9. The left third is just a stack
# of identity matrices ... the function calls a helper function to get the remaining
# two thirds

calc.C <- function(ns){
  
  pt.no = dim(ns)[1]
  C = cbind(kronecker(rep(1,pt.no), diag(c(1,1,1))), C.helper(ns) )
  return(C) 
}

C.helper <- function(ns){
  
  # The first case allows for n to be just one Bloch vector
  if (class(ns) == "numeric"){
    if (length(ns) != 3){
      return(FALSE)
    }
    return(C.helper.helper(ns))    
  }
  else{
    # Check to see ns is a stack of 3-vectors
    if (dim(ns)[2] != 3){
      return(FALSE)
    }
    else{
      mtrx = matrix(0, dim(ns)[1]*3, 6)
      
      # For each Bloch vector, compute the relevant 3x6 matrix
      for (j in 1:dim(ns)[1]){
        mtrx[(j-1)*3+1:3,] = C.helper.helper(ns[j,]) 
      }
      return(mtrx)  
    }
  }
}

# This function calculates the 3x6 matrix piece for each Bloch vector
C.helper.helper <- function(n){
  # This is the middle third 
  piece1 = diag(n) - n %*% t(c(1,1,1))
  # This is the right third
  piece2 = matrix(sum(n),3,3)- n %*% t(c(1,1,1)) - c(1,1,1) %*% t(n) - 
    diag(sum(n),3) + 2*diag(n)
  return(cbind(piece1,piece2)) 
}