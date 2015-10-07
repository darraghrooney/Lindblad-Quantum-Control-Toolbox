# This function generates an arbitrary number of random Lindblad systems, and calcs
# their horizons and optimizing hamiltonians. The rot matrix is ignored for
# simplicity.

data_generate <- function(n_rows){
  
  # Initialize data frame
  test_data = data.frame(a11 = double(n_rows), a22 = double(n_rows),  
          a33 = double(n_rows), a23i = double(n_rows), a31i = double(n_rows),  
          a12i = double(n_rows), nM1 = double(n_rows),nM2 = double(n_rows), 
          nM3 = double(n_rows), hM1 = double(n_rows),hM2 = double(n_rows), 
          hM3 = double(n_rows), zero_check = double(n_rows), 
          CP_check = double(n_rows), hM_check = double(n_rows))

  # Sweep over data points
  for (j in 1:n_rows){

    # Generate random system
    rA = rand_A(FALSE)
    
    # Vectorize the matrix
    test_data[j, 1:6] = A_vec(rA$A)
    
    # Treat special cases
    if (min(rA$a) == max(rA$a)){
        nM = rA$b/2/max(rA$a)
      }
    else if (min(rA$a) == 0){ 
        nM = rA$b/ sum(rA$a) 
    }
    
    # Calculate horizon
    else{
        nM = horizon_find(rA$a, rA$b)
      }
    test_data[j, 7:9] = nM

    # For zero horizon, Hamiltonian can be zero
    if (max(abs(nM)) < 1e-10 ){hM = c(0,0,0)}
    
    # Calculate Hamiltonian
    else {
      hM = -h_mat(nM) %*% ( rA$b + diag(as.vector(rA$a)) %*% nM)/(sum(nM^2))
    }
    test_data[j, 10:12] = hM
      
    # Check that horizon is a zero and a CP
    zero_check = t(nM) %*% ( rA$b + (diag(as.vector(rA$a)) - diag(rep(sum(rA$a)),3))%*% nM )
    deriv = rA$b + 2*diag(as.vector(rA$a)) %*% nM
    if (max(abs(deriv)) < 1e-10 ){CP_check = 0}
    else {
      CP_check = sqrt(sum((h_mat(nM) %*% deriv)^2) / sum(deriv^2))
    }
      
    # Check that Hamiltonian renders horizon stationary
    nderiv = h_mat(hM)%*%nM + rA$b + (diag(as.vector(rA$a)) - diag(rep(sum(rA$a)),3))%*% nM
    hM_check = sum(nderiv^2)
      
    test_data[j,13:15] = c(zero_check, CP_check, hM_check)
  }
  
  return(test_data)
}

# Helper function to vectorize the system matrix
A_vec <- function(A){
  
  Av = rep(0,6)
  Av[1:3] = diag(Re(A))
  Av[4] = Im(A[2,3])
  Av[5] = Im(A[3,1])
  Av[6] = Im(A[1,2])
  return(Av)
  
}

# Helper function to matrixify Hamiltonian
h_mat <- function(hv){
  hm = diag(rep(0,3))
  hm[1,2] = hv[3] 
  hm[2,3] = hv[1] 
  hm[3,1] = hv[2] 
  return( - hm + t(hm) )
}