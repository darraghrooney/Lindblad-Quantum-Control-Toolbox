# This function takes a Lindblad system and a Hamiltonian and assess how close
# the Hamiltonian is to the system's maximizing Hamiltonian. There are six 
# measures: relative error in the max Ham, relative error in the horizon,
# relative error in the purity, absolute error in the purity, absolute
# purity gained compared to the zero Hamiltonian, and the relative purity gained.

# is.rot indicates whether Re(A) is undiagonalized. Diagonalized Re(A) makes
# it faster to process
hM_check <- function(A, hM_to_check, is.rot=FALSE){

  # Breakdown system if undiagonalized Re(A)
  if (is.rot){
    A_full = A_breakdown(A)
    rot = A_full$rot
    a = A_full$a
    b = A_full$b
  }
  else{
    a = sort(diag(Re(A)), decreasing=TRUE)
    b = b_vec(Im(A))
    rot = diag(c(1,1,1))    
  }

  # Calculate horizon and maximizing Hamiltonian 
  nM = horizon_find(a,b)
  hM = -h_mat(nM) %*% ( b + diag(a) %*% nM)/(sum(nM^2))

  # Calculate horizons from input and zero Hamiltonians 
  nM_to_check = solve(h_mat(hM_to_check) + Re(A) - sum(diag(Re(A)))*diag(c(1,1,1)), -b)
  n_zero = solve(Re(A) - sum(diag(Re(A)))*diag(c(1,1,1)), -b)
  
  # Calculate errors
  hM_error = sqrt(sum(hM - hM_to_check)^2)/sqrt(sum(hM^2))
  nM_error = sqrt(sum(nM - nM_to_check)^2)/sqrt(sum(nM^2))
  lost_purity_abs = sqrt(sum(nM^2)) - sqrt(sum(nM_to_check^2))
  lost_purity_rel = lost_purity_abs/sqrt(sum(nM^2))
  purity_gain_abs = sqrt(sum(nM_to_check^2)) - sqrt(sum(n_zero^2))
  purity_gain_rel = purity_gain_abs/(sqrt(sum(nM^2)) 
                                     - sqrt(sum(n_zero^2)))

  return(list(hM_error,nM_error,lost_purity_rel,lost_purity_abs, 
              purity_gain_abs,purity_gain_rel))
}

# Helper function to matrixify Hamiltonian
h_mat <- function(hv){
  hv = as.numeric(hv)
  hm = diag(rep(0,3))
  hm[1,2] = hv[3] 
  hm[2,3] = hv[1] 
  hm[3,1] = hv[2] 
  return( - hm + t(hm) )
}

# Helper function to vectorize the imag part of the A matrix
b_vec <- function(bm){
  bv = as.matrix(rep(0,3))
  bv[1] = 2*bm[2,3]
  bv[2] = 2*bm[3,1]
  bv[3] = 2*bm[1,2]
  return( bv )
}