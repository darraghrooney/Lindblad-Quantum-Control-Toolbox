# This function generates a random n=2 Lindblad system
# By default the a-vector and b-vector are generated randomly.
# If the rot_switch is set to true, the system is given a random
# SO(3) rotation

rand_A <- function( rot_switch = FALSE){

  # a has three uniform rvs on [0,100], and is sorted 
  a = matrix(sort(c(100*runif(3)), decreasing=TRUE))

  # To get b, we pick a uniform variable in the unit ball, and stretch
  # accordingly
  cosph = 2*runif(1)-1
  theta = 2*pi*runif(1)
  brad = runif(1)^(1/3)
  
  bangles = c(cos(theta)*sin(acos(cosph)), sin(theta)*sin(acos(cosph)), cosph)
  bax = c(sqrt(a[2]*a[3]),sqrt(a[3]*a[1]),sqrt(a[1]*a[2]))
  b = 2*brad*bax*bangles
  
  # Calculate the a matrix
  A = diag(as.vector(a)) + 1i*b_mat(b)
  
  # Calculate the rotation matrix if desired. Uses a method 
  # involving quaternions
  rot = diag(rep(1,3))
  if(rot_switch){
    u = runif(3)
    q = c(sin(2*pi*u[2]), cos(2*pi*u[2]),sin(2*pi*u[3]),cos(2*pi*u[3]))
    q[1:2] = q[1:2]*sqrt(1-u[1])
    q[3:4] = q[3:4]*sqrt(u[1])
    rot[1,1] = rot[1,1] - 2*q[3]^2 - 2*q[4]^2
    rot[2,2] = rot[2,2] - 2*q[4]^2 - 2*q[2]^2
    rot[3,3] = rot[3,3] - 2*q[2]^2 - 2*q[3]^2
    rot[1,2] = 2*q[2]*q[3] - 2*q[1]*q[4]
    rot[2,1] = 2*q[2]*q[3] + 2*q[1]*q[4]
    rot[2,3] = 2*q[3]*q[4] - 2*q[1]*q[2]
    rot[3,2] = 2*q[3]*q[4] + 2*q[1]*q[2]
    rot[3,1] = 2*q[2]*q[4] - 2*q[1]*q[3]
    rot[1,3] = 2*q[2]*q[4] + 2*q[1]*q[3]

  }
  
  return(list(A = rot %*% A %*% t(rot), a=a, b=b, rot=rot))
}

# Helper function that converts a b-vector into its matrix form
b_mat <- function(bv){
  bm = diag(rep(0,3))
  bm[1,2] = bv[3]/2 
  bm[2,3] = bv[1]/2 
  bm[3,1] = bv[2]/2 
  return( bm - t(bm) )
}

# Inverse of the b_mat function
b_vec <- function(bm){
  bv = as.matrix(rep(0,3))
  bv[1] = 2*bm[2,3]
  bv[2] = 2*bm[3,1]
  bv[3] = 2*bm[1,2]
  return( bv )
}