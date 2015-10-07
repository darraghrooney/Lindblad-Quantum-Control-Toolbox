# This function generates the horizons and necessary Hamiltonians by sweeping
# over the parameter space. The parameter space is nine-d, but only 4 directions
# are necessary. Three can be rotated away (this function assumes this has been 
# done.) The largest a is taken to be 100, since multiplying all parameters by 
# a constant does not affect the data. Additionally, the magnitude of b_scal is 
# taken to be one, since the Hamiltonian is invariant under different magnitudes.

# Input variables: the discretization parameters (maximum indices) for a- 
# and b- space. Also the starting and final indices for a2 and a3 are entered so 
# that one canselect portions of the parameter space.

landscape_generate <- function(a_disc, b_disc, aj2i, aj2f, aj3i, aj3f){
  
  # Initialization  
  no_b = (b_disc+1)^2
  tri_dim = max(0,aj3f - aj2i+1)
  n_rows = no_b*( (aj2f-aj2i+1)*(aj3f-aj3i+1) - (tri_dim)*(tri_dim-1)/2  ) 
  sys_data = data.frame(a2_ind = double(n_rows),          
      a3_ind = double(n_rows), cosphi_ind = double(n_rows), theta_ind = 
      double(n_rows), nM1 = double(n_rows),nM2 = double(n_rows), nM3 = 
      double(n_rows), hM1 = double(n_rows),hM2 = double(n_rows), hM3 = 
      double(n_rows), zero_check = double(n_rows), CP_check = double(n_rows), 
      hM_check = double(n_rows))
  count = 0
  a1 = 100;
  
  # Sweep over a-space
  for (j1 in aj2i:aj2f){
  for (j2 in aj3i:min(j1,aj3f,a_disc)){
    a2 = 100*j1/a_disc
    a3 = 100*j2/a_disc
    a = c(a1,a2,a3)
    
    # Sweep over b-space (actually an eighth of b-space is necessary due to symmetry)
    for (k1 in 0:b_disc){
    for (k2 in 0:b_disc){
      
      # Calculate b      
      sk1 = sqrt(1-(k1/b_disc)^2)
      b1 = 2*sqrt(a2*a3)*cos(pi*k2/b_disc/2)*sk1
      b2 = 2*sqrt(a3*a1)*sin(pi*k2/b_disc/2)*sk1
      b3 = 2*sqrt(a1*a2)*k1/b_disc
      b = matrix(c(b1,b2,b3))

      # Calculate horizon
      if(j2==a_disc){nM = b/200}
      else if(j2==0) { nM = b/(a1+a2) }
      else{nM = horizon_find(a,b)}

      # Calculate Hamiltonian
      if (max(abs(nM)) < 1e-10 ){hM = c(0,0,0)}
      else {
        hM = -h_mat(nM) %*% ( b + diag(a) %*% nM)/(sum(nM^2))
      }
      
      # Check that horizon is a zero and a CP
      zero_check = t(nM) %*% ( b + (diag(a) - diag(rep(sum(a)),3))%*% nM )
      deriv = b + 2*diag(a) %*% nM
      if (max(abs(deriv)) < 1e-10 ){CP_check = 0}
      else {
        CP_check = sqrt(sum((h_mat(nM) %*% deriv)^2) / sum(deriv^2))
      }
      
      # Check that Hamiltonian renders horizon stationary
      nderiv = h_mat(hM)%*%nM + b + (diag(a) - diag(rep(sum(a)),3))%*% nM
      hM_check = sum(nderiv^2)
      count = count + 1
      
      # Add data to data frame
      sys_data[count,] = c(j1,j2,k1,k2,nM,hM,zero_check,CP_check,hM_check)
    }}}}
  
  
  return(sys_data)
}

# Helper function to matrixify Hamiltonian
h_mat <- function(hv){
  hm = diag(rep(0,3))
  hm[1,2] = hv[3] 
  hm[2,3] = hv[1] 
  hm[3,1] = hv[2] 
  return( - hm + t(hm) )
}