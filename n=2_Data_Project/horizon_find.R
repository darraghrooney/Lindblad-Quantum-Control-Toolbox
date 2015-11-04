# This function calculates the optimal horizon for an n=2 dissipative quantum system.
# Input is a 3-vector a, which is non-negative, and a 3-vector b, which obeys the 
# inequality a1 b1^2 + a2 b2^2 + a3 b3^3 <= 4 a1 a2 a3. We also assume a1 >= a2 >= a3.

# Output are five numbers. The first is the Lagrange multiplier used to find the
# horizon, the next three are the co-ordinates of the horizon, and the fifth is the
# horizon magnitude.

horizon_find <- function(a, b){
  
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
      mus = max(a)*find_mus3(a[c(J1,J2,J3)]/max(a),b[c(J1,J2,J3)]/max(a))
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
    mus = max(a)*find_mus5(a/max(a),b/max(a))
    
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

# Helper function to matrixify vector for cross product
cross_mat <- function(hv){
  hm = diag(rep(0,3))
  hm[1,2] = hv[3] 
  hm[2,3] = hv[1] 
  hm[3,1] = hv[2] 
  return( - hm + t(hm) )
}

# Calculates the Lagrangian multiplier for a Lindblad system
# with two non-zero bj's. Involves solving a cubic equation.
# It is expected that A1 is on the order of 100 ... coefficients
# are re-scaled to maintain accuracy

find_mus3 <- function(a,b){
  
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

find_mus5 <- function(a,b){

  # Calculate coefficients, including re-scaling factors
  coeff_exp = matrix(0,6,3);
  for(j in 1:3){
    j1 = (j %% 3) + 1;
    j2 = ((j+1) %% 3)+1;
    A1 = a[j]; B1 = b[j];
    A2 = a[j1]; B2 = b[j1];
    A3 = a[j2]; B3 = b[j2];

    coeff_exp[1,j] = 2;   
    coeff_exp[2,j] = -4*(A2+A3);   
    coeff_exp[3,j] = 2*(A2^2+A3^2+4*A2*A3);   
    coeff_exp[4,j] = -4*(A2*A3^2+A3*A2^2);   
    coeff_exp[5,j] = 2*A2^2*A3^2;   
    coeff_exp[2:6,j] = coeff_exp[2:6,j] - coeff_exp[1:5,j]*(A2/2+A3/2+A1);
    coeff_exp[,j] = coeff_exp[,j]*B1^2;
  }
  coeff = rowSums(coeff_exp)
  
  # Solve polynomial and scale back
  mus = sum(a) - polyroot(coeff[6:1])

  # Before returning, toss out non-real solutions and sort
  # in descending order
  return( sort(Re(mus[abs(Im(mus)) < 1e-3]), decreasing=TRUE) )
}

# Helper function that calculates the components of b
# from its angle parametrization
b_from_angles <- function(a,cosphi,theta,B){
  
  b=rep(0,3)
  b[1] = 2*sqrt(a[2]*a[3])*B*cos(theta)*sin(acos(cosphi))
  b[2] = 2*sqrt(a[1]*a[3])*B*sin(theta)*sin(acos(cosphi))
  b[3] = 2*sqrt(a[2]*a[1])*B*cosphi
  
  return(b)
}

# Given a point in a Bloch ball for a Lindblad system, computes
# the necessary (non-unique) Hamiltonian
h_from_n <- function(n,a,b){
  
  # Matrixify n
  n_mat = cross_mat(n)
  
  # Compute Hamiltonian
  if (sum(n^2) < 1e-8){
    return(rep(0,3))
  }
  hM = - n_mat %*% ( b + diag(a) %*% n)/(sum(n^2))
  return(hM)
}