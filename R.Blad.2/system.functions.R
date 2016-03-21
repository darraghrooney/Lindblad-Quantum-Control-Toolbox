# List of functions in this file:
#   rand.A()
#   b.from.angles()
#   rot.from.u()
#   u.from.rot()
#   cross.mat()
#   cross.vec()
#   A.breakdown()
#   A.assemble()
#   h.from.n()
#   n.from.h()
#   horizon.find()
#   find.mus5()
#   find.mus3()
#   mag()
#   tr()
#   perp.basis()

##############################################################

# This function generates a random n=2 Lindblad system.
# By default the a-vector and b-vector are generated randomly.
# If the rot.switch is set to true, the system is given a random
# SO(3) rotation

rand.A <- function( rot.switch = FALSE ){
  
  # a is a sorted vector whose largest element is 1. The other two are uniform on [0,1]
  a = c(1,sort(runif(2),decreasing=TRUE))
  
  # To get b, we pick a uniform point in the 1/8th unit ball, and stretch accordingly
  cosphi = runif(1)
  theta = pi*runif(1)/2
  brad = runif(1)^(1/3)
  b = b.from.angles(a,cosphi,theta,brad)
  
  # Calculate the A matrix
  A = diag(as.vector(a)) - 1i*cross.mat(b)/2
  
  # Calculate the rotation matrix if desired. Uses a method involving quaternions
  rot = diag(rep(1,3))
  u = c(0,.25,0)
  if(rot.switch){
    u = runif(3)
    rot = rot.from.u(u)
    u = u.from.rot(rot)   # To eliminate ambiguities, reset u to the given output from u.from.rot
  }
  
  return(list(A = rot %*% A %*% t(rot), a=a, b=b, u=u))
}

# Helper function that calculates the components of b from its angle parametrization
b.from.angles <- function(a,cosphi,theta,B){
  
  b=rep(0,3)
  b[1] = 2*sqrt(a[2]*a[3])*B*cos(theta)*sin(acos(cosphi))
  b[2] = 2*sqrt(a[1]*a[3])*B*sin(theta)*sin(acos(cosphi))
  b[3] = 2*sqrt(a[2]*a[1])*B*cosphi
  
  return(b)
}

# Helper function that maps a triple on [0,1] to a rotation in SO(3)
rot.from.u <- function(u){

  rot = diag(rep(1,3))
  
  # Map input to a quaternion
  q = c(sin(2*pi*u[2]), cos(2*pi*u[2]),sin(2*pi*u[3]),cos(2*pi*u[3]))
  q[1:2] = q[1:2]*sqrt(1-u[1])
  q[3:4] = q[3:4]*sqrt(u[1])

  # Map quaternion to rotation
  rot[1,1] = 1 - 2*q[3]^2 - 2*q[4]^2
  rot[2,2] = 1 - 2*q[4]^2 - 2*q[2]^2
  rot[3,3] = 1 - 2*q[2]^2 - 2*q[3]^2
  rot[1,2] = 2*q[2]*q[3] - 2*q[1]*q[4]
  rot[2,1] = 2*q[2]*q[3] + 2*q[1]*q[4]
  rot[2,3] = 2*q[3]*q[4] - 2*q[1]*q[2]
  rot[3,2] = 2*q[3]*q[4] + 2*q[1]*q[2]
  rot[3,1] = 2*q[2]*q[4] - 2*q[1]*q[3]
  rot[1,3] = 2*q[2]*q[4] + 2*q[1]*q[3]
  
  return(rot)  
}

# Inverse of rot.from.u. Since quarternion->rotation is not injective, answer is not unique.

u.from.rot <- function(rot){
  
  u = rep(0,3)
  q = rep(0,4)
  
  # General case has trace(rot) not equal to -1
  if ( abs(tr(rot) + 1) > 1e-10) {
    
    # Compute quaternion
    q[1] = sqrt( 1 + tr(rot) )/2      # Can be +/-, but not zero. Choose +.
    q[2] = (rot[3,2]-rot[2,3])/4/q[1]
    q[3] = (rot[1,3]-rot[3,1])/4/q[1]
    q[4] = (rot[2,1]-rot[1,2])/4/q[1]
  
    # Compute u from quaternion
    u[1] = q[3]^2+q[4]^2  
    if ( u[1] != 1){
      u[2] = atan2(q[1], q[2])/2/pi   # If u[1] is one, u[2] is free. Choose zero. 
    }
    if ( u[1] != 0){
      u[3] = atan2(q[3], q[4])/2/pi   # If u[1] is zero, u[3] is free. Choose zero.
    }
  }
  
  # If trace(rot) is -1, and all off-diags are non-zero.
  else if ( min(abs( c(rot[1,2],rot[2,3],rot[3,1]) )) > 1e-14)
  {
    # q[1] = 0. Compute remaining q's. There is a double sign ambiguity. 
    q[2] = sqrt(rot[1,2]*rot[1,3]/rot[2,3]/2)*sign(rot[2,3])
    q[3] = sqrt(rot[1,2]*rot[2,3]/rot[1,3]/2)*sign(rot[3,1])
    q[4] = sqrt(rot[2,3]*rot[1,3]/rot[1,2]/2)*sign(rot[1,2])
    
    # Compute u from quaternion. No ambiguity.
    u[1] = q[3]^2+q[4]^2
    u[2] = atan2(q[1],q[2])/2/pi
    u[3] = atan2(q[3],q[4])/2/pi
  }
  
  # If trace(rot) is -1, and all off-diags are zero, diag has one 1, and two -1's.
  # If rot[1,1] = 1, u[1] = 0, u[2] = 0 or 1, u[3] is free. Choose all zeros, so do nothing.
  
  else if ( max(abs( c(rot[1,2],rot[2,3],rot[3,1]) )) < 1e-14)
  {
    # Case where rot[2,2] = 1. u[2] is free, choose zero. u[3] is 1/4 or 3/4. Choose 1/4.
    if ( rot[2,2] > 0){ u[1] = 1; u[3] = 1/4 }
    # Case where rot[3,3] = 1. u[2] is free, choose zero. u[3] is 0, 1/2 or 1. Choose 0.
    else if (rot[3,3] > 0){ u[1] = 1 }
  }

  # If trace(rot) is -1 and the off-diags are partially but not completely zero, then
  # exactly four out of the six are zero, and rot is symmetric.
  else if ( abs(rot[1,2]) > 1e-14){
    q[2] = sqrt((1 - rot[2,2])/2)               # Sign is ambiguous. Choose positive.
    q[3] = sqrt((1 - rot[1,1])/2)*sign(rot[1,2])
    u = c( q[3]^2+q[4]^2, atan2(q[1],q[2])/2/pi, atan2(q[3],q[4])/2/pi)
  }
  else if ( abs(rot[2,3]) > 1e-14){
    q[3] = sqrt((1 - rot[3,3])/2)         # Sign is ambiguous. Choose positive.
    q[4] = sqrt((1 - rot[2,2])/2) * sign(rot[2,3])
    u = c(1, 0, atan2(q[3],q[4])/2/pi)    # u[2] is free. Choose zero.
  }
  else if ( abs(rot[3,1]) > 1e-14){
    q[4] = sqrt((1 - rot[1,1])/2)         # Sign is ambiguous. Choose positive.
    q[2] = sqrt((1 - rot[3,3])/2) * sign(rot[3,1])
    u = c( q[3]^2+q[4]^2, atan2(q[1],q[2])/2/pi, atan2(q[3],q[4])/2/pi)    
  }
  
  # The atan2 returns values in (-pi,pi], and we want [0,2pi), so adjust.
  u[2:3] = u[2:3] - floor(u[2:3])
  return(u)
}

# Helper function to matrixify vector for cross product
cross.mat <- function(vc){
  cmt = diag(rep(0,3))
  cmt[1,2] = vc[3] 
  cmt[2,3] = vc[1] 
  cmt[3,1] = vc[2] 
  return( - cmt + t(cmt) )
}

# Inverse of the cross.mat function
cross.vec <- function(cmt){
  cv = as.matrix(rep(0,3))
  cv[1] = cmt[3,2]
  cv[2] = cmt[1,3]
  cv[3] = cmt[2,1]
  return( cv )
}

# This function takes a positive-definite 3x3 complex matrix, and decomposes it
# into three parts: the sorted vector of eigenvalues of Re(A), the eigenvector
# matrix of Re(A), and the b-vector representing the imaginary parts of the
# off-diagonal elements. The rotation is chosen so that the components of b are 
# all non-negative. If b has zero components, there is a binary ambiguity in the
# relative axis, but there's no systematic way to resolve it.

A.breakdown <- function(A){
  
  eiginfo = eigen(Re(A), symmetric=TRUE)
  a = eiginfo$values
  rot = eiginfo$vectors*det(eiginfo$vectors)
  b = 2 * t(rot) %*% cross.vec(-Im(A))
  rot = rot %*% diag( as.vector(sign(b) > -.5)*2 -1 )  # Flip column vector so that b has non-negative components
  u = u.from.rot(rot)
  
  return(list(a=a,b=abs(b), u=u ))
}

# Inverse function of A.breakdown(). If b has zero components, A.breakdown(A.assemble())
# may not be the identity.

A.assemble <- function(Ab){
  rot = rot.from.u(Ab$u)
  return( rot %*% ( diag(as.vector(Ab$a)) - 1i/2*cross.mat(Ab$b)) %*% t(rot) )
}

# Given a point in a Bloch ball for a Lindblad system, computes
# the necessary (non-unique) Hamiltonian
h.from.n <- function(n,a,b){
  
  # Compute Hamiltonian
  if (sum(n^2) < 1e-8){
    return(rep(0,3))
  }
  hM = - cross.mat(n) %*% ( b + diag(a) %*% n)/(mag(n)^2)
  return(hM)
}

# Given a Hamiltonian and Lindblad system, return the asymptotic n.
# If the system is degenerate, return FALSE.
n.from.h <- function(h,a,b){
  
  # Compute dynamic matrix
  fullmat = cross.mat(h) + diag(a - sum(a))
  
  n = try(solve(fullmat, -b), TRUE)
  if (class(n) == "try-error"){
    n = c(0,0,0)
  }
  return(n)
}

# This function calculates the optimal horizon for an n=2 dissipative quantum system.
# Input is a 3-vector a, which is non-negative, and a 3-vector b, which obeys the 
# inequality a1 b1^2 + a2 b2^2 + a3 b3^3 <= 4 a1 a2 a3. We also assume a1 >= a2 >= a3.

# Output are five numbers. The first is the Lagrange multiplier used to find the
# horizon, the next three are the co-ordinates of the horizon, and the fifth is the
# horizon magnitude.

horizon.find <- function(a, b){
  
  horizon = rep(0,5) 
  
  # If b has only one non-zero component, we know the analytic solution 
  if (sort(abs(b))[2] < 1e-10){
    J1 = which.max(abs(b))        # Find non-zero component
    horizon[1] = (sum(a) - a[J1])/2         # Solution
    horizon[1+J1] = b[J1]/(sum(a)-a[J1])  # Solution
    
    J2 = setdiff(1:3,c(1,J1))
    if (J1 != 1 && a[1] > 2*a[J1]+a[J2]){   # There is an alternate solution
      horizon[1] = sum(a) - a[1]              # which may beat the standard solution
      horizon[1+J1] = b[J1]/2/(a[1]-a[J1])
      horizon[2] = horizon[1+J1]*sqrt((a[1]-2*a[J1]-a[J2])/(sum(a)-a[1]))
    }
    horizon[5] = sqrt(sum(horizon[2:4]^2))  # Calculate magnitude
  }
  
  # If b has two non-zero components, we also have a special case.
  else if (sort(abs(b))[1] < 1e-10){
    J1 = which.min(abs(b))
    J2 = (J1 %% 3) + 1
    J3 = ((J1+1) %% 3)+1
    
    # The standard solution involves solving a cubic polynomial which we 
    # assign to a helper function. The data is re-scaled so that a1=1. This
    # avoids very large polynomial co-efficients.
    mus = max(a)*find.mus3(a[c(J1,J2,J3)]/max(a),b[c(J1,J2,J3)]/max(a))
    for (k in 1:length(mus)){
      nM = b/2/(-mus[k]+sum(a)-a)
      if (min(is.finite(nM))==0){next}
      nMmag2 = sum(nM^2)
      if (nMmag2 > horizon[5]^2 && nMmag2 - 1 < 1e-10){
        horizon = c(mus[k],nM,sqrt(nMmag2))
      }
    }
    # We must construct the possible alternate solution
    nMsp = rep(0,3)
    nMsp[J2] = b[J2]/2/(a[J1]-a[J2])
    nMsp[J3] = b[J3]/2/(a[J1]-a[J3])
    spv = nMsp[J2]^2 *(a[J1]-2*a[J2]-a[J3])/(sum(a)-a[J1]) + 
      nMsp[J3]^2 *(a[J1]-2*a[J3]-a[J2])/(sum(a)-a[J1])
    
    # Check to see if alternate solution is valid and optimal
    if (spv > 0 && spv-1 < 1e-10){  # Valid?
      nMsp[J1] = sqrt(spv)
      nMspmag = sqrt(sum(nMsp^2)) 
      if (nMspmag > horizon[5]){    # Optimal?
        horizon = c(sum(a) - a[J1], nMsp, nMspmag) # Update if appropriate
      }
    }
  }
  # In general, we must solve a quintic polynomial
  else{
    # Use helper function to solve polynomial. Rescale to avoid numeric craziness
    mus = max(a)*find.mus5(a/max(a),b/max(a))
    
    # See which solution is best
    for (k in 1:length(mus)){
      nM = b/2/(-mus[k]+sum(a)-a)
      if (min(is.finite(nM))==0){next}
      nMmag2 = sum(nM^2)
      if (nMmag2 > horizon[5]^2 && nMmag2 - 1 < 1e-10){   # Optimal and valid?
        horizon = c(mus[k],nM,sqrt(nMmag2))
      }
    }
  }
  return(horizon)
}

# Calculates the Lagrangian multiplier for a Lindblad system
# with two non-zero bj's. Involves solving a cubic equation.
# It is expected that A1 is on the order of 100 ... coefficients
# are re-scaled to maintain accuracy

find.mus3 <- function(a,b){
  
  A1 = a[1]
  A2 = a[2]
  A3 = a[3]
  B2 = b[2]
  B3 = b[3]
  
  # Calculate polynomial co-efficients. There are scaling factors in there
  
  coeff = c( -2*(B2^2+B3^2),
             B2^2*(4*A2+5*A1+A3)+B3^2*(4*A3+5*A1+A2),
             -2*B2^2*(A1+A2)*(2*A1+A2+A3)-2*B3^2*(A1+A3)*(2*A1+A2+A3),
             B2^2*(A1+A3)*(A1+A2)^2+B3^2*(A1+A2)*(A1+A3)^2)
  # Solve polynomial and scale back  
  mus = polyroot(coeff[4:1])
  
  # Before returning, toss out non-real solutions and sort
  # in descending order
  return( sort(Re(mus[abs(Im(mus)) < 1e-4]), decreasing=TRUE) )
}

# Calculates the Lagrangian multiplier for a Lindlbad system
# with three non-zero bj's. Involves solving a quintic equation.
# It is expected that max(a) is on the order of 100, as polynomial
# co-efficients are re-scaled to maintain accuracy

find.mus5 <- function(a,b){
  
  # Calculate coefficients, including re-scaling factors
  coeff.exp = matrix(0,6,3);
  for(j in 1:3){
    j1 = (j %% 3) + 1;
    j2 = ((j+1) %% 3)+1;
    A1 = a[j]; B1 = b[j];
    A2 = a[j1]; B2 = b[j1];
    A3 = a[j2]; B3 = b[j2];
    
    coeff.exp[1,j] = 2;   
    coeff.exp[2,j] = -4*(A2+A3);   
    coeff.exp[3,j] = 2*(A2^2+A3^2+4*A2*A3);   
    coeff.exp[4,j] = -4*(A2*A3^2+A3*A2^2);   
    coeff.exp[5,j] = 2*A2^2*A3^2;   
    coeff.exp[2:6,j] = coeff.exp[2:6,j] - coeff.exp[1:5,j]*(A2/2+A3/2+A1);
    coeff.exp[,j] = coeff.exp[,j]*B1^2;
  }
  coeff = rowSums(coeff.exp)
  
  # Solve polynomial and scale back
  mus = sum(a) - polyroot(coeff[6:1])
  
  # Before returning, toss out non-real solutions and sort
  # in descending order
  return( sort(Re(mus[abs(Im(mus)) < 1e-3]), decreasing=TRUE) )
}

# Apparently, R has no function for computing the magnitude of a vector (!?), 
# so here it is.

mag <- function(n){
  if (class(n) == "data.frame"){
    return(sqrt(rowSums(n^2)))
  }
  else{
    return(sqrt(sum(n^2)))
  }
}

# R also has no matrix trace function.

tr <- function(A){
  return(sum(diag(A)))
}

# This is a function that, given a vector in R3, returns an orthonormal basis for
# the orthogonal space.

perp.basis <- function(v){
  
  v = v/mag(v)
  w1 = cross.mat(t(c(1,0,0))) %*% v
  w2 = cross.mat(t(c(0,1,0))) %*% v
  w3 = cross.mat(t(c(0,0,1))) %*% v
  evs = eigen(cbind(w1,w2,w3))$vectors[,1]*sqrt(2)
  return(cbind(Re(evs), Im(evs)))
}