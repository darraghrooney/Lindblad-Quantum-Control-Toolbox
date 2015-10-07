# This function takes an arbitrary Lindblad system, and uses the landscape
# to determine the system's maximizing Hamiltonian, using interpolation. 

landscape_interpolate <- function(A, need_to_read){
 
  # Read in landscape if necessary
  if (need_to_read){
    library(R.matlab)
    rawls = readMat("landscape.MAT")
    landscape = data.frame(rawls$landscape[,,1])
  }
  
  # Breakdown system matrix into relevant parts
  A_full = A_breakdown(A)
  a = A_full$a
  b = A_full$b
  rot = A_full$rot
  
  # Normalize system parameters
  b_sgn = sign(b)
  anorm = 100/a[1]*a
  bnorm = 100/a[1]*abs(b)
  
  # Extract b-angles from system parameters
  if (anorm[2] < 1e-10){
    bn_cosphi=0
    bn_theta=0}
  else if (anorm[3] < 1e-10){
    bn_theta = 0
    bn_cosphi = sign(bnorm[3])
  }
  else{
    bn_cosphi = bnorm[3]/2/sqrt(anorm[1]*anorm[2])
    bn_theta = atan2(bnorm[2]/2/sqrt(anorm[1]*anorm[3]), 
                     bnorm[1]/2/sqrt(anorm[2]*anorm[3]))
  }
    
  # Fetch the necessary variables in the landscape
  land_a = data.frame(a2_ind=landscape$a2_ind, a3_ind=landscape$a3_ind)
  land_ang = data.frame(cosphi_ind=landscape$cosphi_ind, 
                        theta_ind=landscape$theta_ind)
  land_hM = data.frame(hM1=landscape$hM1,hM2=landscape$hM2, hM3=landscape$hM3)
  
  # Get some parameters from landscape
  land_a_samp = max(land_a$a2_ind)
  del_a = 100/land_a_samp
  land_b_samp = max(land_ang$cosphi_ind)
  
  # Determine which grid points should be used for interpolation
  fract_ind = c(anorm[2]/del_a ,anorm[3]/del_a,
                bn_cosphi*land_b_samp, bn_theta*land_b_samp*2/pi)
  root_ind = round( fract_ind)
  if (root_ind[1]==root_ind[2]){
    root_ind[1] = ceiling(anorm[2]/del_a ) 
    root_ind[2] = floor(anorm[3]/del_a )
  }
  ind_signs = sign(fract_ind - root_ind)

  # Calculate the Hamiltonians for the relevant grid points
  h_root = landscape[ landscape$a2_ind == root_ind[1] & 
    landscape$a3_ind == root_ind[2] & landscape$cosphi_ind == root_ind[3] & 
    landscape$theta_ind == root_ind[4], c(8:10)]
  h_a2 = landscape[ landscape$a2_ind == root_ind[1] + ind_signs[1] & 
    landscape$a3_ind == root_ind[2] & landscape$cosphi_ind == root_ind[3] & 
    landscape$theta_ind == root_ind[4], c(8:10)]
  h_a3 = landscape[ landscape$a2_ind == root_ind[1] & landscape$a3_ind == 
    root_ind[2]  + ind_signs[2] & landscape$cosphi_ind == root_ind[3] & 
    landscape$theta_ind == root_ind[4], c(8:10)]
  h_cp = landscape[ landscape$a2_ind == root_ind[1] & 
    landscape$a3_ind == root_ind[2] & landscape$cosphi_ind == root_ind[3] 
    + ind_signs[3] & landscape$theta_ind == root_ind[4], c(8:10)]
  h_th = landscape[ landscape$a2_ind == root_ind[1] & 
    landscape$a3_ind == root_ind[2] & landscape$cosphi_ind == root_ind[3] & 
    landscape$theta_ind == root_ind[4] + ind_signs[4], c(8:10)]
  
  # Interpolate between calculate Hamiltonians
  hM_interp = h_root + (h_a2-h_root)*abs(fract_ind[1]-root_ind[1])  + 
    (h_a3-h_root)*abs(fract_ind[2]-root_ind[2]) + abs(h_cp-h_root)*(fract_ind[3]- 
                  root_ind[3]) + abs(h_th-h_root)*(fract_ind[4]-root_ind[4])

  # Add necessary signs
  h_sgn = c(b_sgn[2]*b_sgn[3], b_sgn[3]*b_sgn[1], b_sgn[1]*b_sgn[2])
  hM_interp = matrix(as.numeric(hM_interp))

  # Un-normalize and rotate if necessary
  return(rot%*% matrix(h_sgn*hM_interp*a[1]/100))
}
  