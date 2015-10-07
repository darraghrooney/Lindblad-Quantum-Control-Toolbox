# This is a function to estimate maximizing Hamiltonians and testing them.
# The test data frame has 20,000 data points, so "data_selection" should be a 
# subset of 1:20000.

interpolation_test <- function(data_selection, need_to_read=FALSE){
  
  # Read in data if necessary
  if (need_to_read){
    library(R.matlab)
    rawls = readMat("landscape.MAT")
    landscape = data.frame(rawls$landscape[,,1])
    rawtd = readMat("test_data.MAT")
    test_data = data.frame(rawtd$test_data[,,1])
  }
  
  # Initialize results data frame
  n_rows = length(data_selection)
  test_results = data.frame(a2norm = double(n_rows), a3norm = double(n_rows),
      cosphi = double(n_rows), theta = double(n_rows),
      hM1 = double(n_rows), hM2 = double(n_rows), hM3 = double(n_rows),
      hM_error = double(n_rows), nM_error = double(n_rows),
      purity_lost_abs = double(n_rows), purity_lost_rel = double(n_rows),
      purity_gained_abs = double(n_rows), purity_gained_rel = double(n_rows)
    )
  
  # Sweep over data points
  count = 1
  for (j in data_selection){

    # Extract data point
    data_to_check = test_data[j,1:15]

    # Re-shape system parameters
    A = A_form(data_to_check[1:6])
    A_break = A_breakdown(A)
    a2norm = A_break$a[2]/A_break$a[1]*100
    a3norm = A_break$a[3]/A_break$a[1]*100
    b1n = A_break$b[1]/2/sqrt(A_break$a[2]*A_break$a[3])
    b2n = A_break$b[2]/2/sqrt(A_break$a[3]*A_break$a[1])
    b3n = A_break$b[3]/2/sqrt(A_break$a[1]*A_break$a[2])
    cosphi = abs(b3n/sqrt(b1n^2+b2n^2+b3n^2))
    theta = acos(abs(cos(atan2(b2n, b3n))))

    # Estimate maximizing Hamiltonian 
    hM_to_check = landscape_interpolate(A, need_to_read)

    # Assess the estimated Hamiltonian
    errors = hM_check(A, hM_to_check)
  
    test_results[count,] = c(a2norm, a3norm, cosphi,theta, hM_to_check, errors[[1]], 
       errors[[2]], errors[[3]], errors[[4]], errors[[5]], errors[[6]])
    count = count + 1
    }
  return(test_results)
}

# Helper function to matrixify the system matrix
A_form <- function(Av){
  
  A= matrix(0,3,3)
  
  A[1,1] = as.numeric(Av[1])
  A[2,2] = as.numeric(Av[2])
  A[3,3] = as.numeric(Av[3])
  A[2,3] = 1i*as.numeric(Av[4])
  A[3,1] = 1i*as.numeric(Av[5])
  A[1,2] = 1i*as.numeric(Av[6])
  A[3,2] = Conj(A[2,3])
  A[1,3] = Conj(A[3,1])
  A[2,1] = Conj(A[1,2])
  
  return(A)
}