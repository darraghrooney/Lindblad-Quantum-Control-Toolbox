# This function takes the parameters for an n=2 Lindblad quantum system
# and calculates the (outermost) horizon in the Bloch ball
# a should be a vector of three non-negative numbers
# b should be a vector of three reals that obey the inequality
# a1 b1^2 + a2 b2^2 + a3 b3^2 <= 4 a1 a2 a3

horizon_find <- function(a,b){
  
  # For zero b, horizon is zero
  if (max(abs(b)) < 1e-10){return( c(0,0,0) )}          

  # For b with only one non-zero element, we know the formula
  if ( max(abs(b[c(1,2)])) < 1e-10 | max(abs(b[c(1,3)])) < 1e-10 |
                max(abs(b[c(2,3)])) < 1e-10 ){
        horizon = c(0,0,0)            
        if (abs(b[1])>1e-10){ horizon[1]=b[1]/(a[2]+a[3])}
        if (abs(b[2])>1e-10){ horizon[2]=b[2]/(a[1]+a[3])}
        if (abs(b[3])>1e-10){ horizon[3]=b[3]/(a[1]+a[2])}
        
        return(horizon)            
      }

  # For b with two non-zero elements, we must solve a cubic polynomial
  if ( min(abs(b)) < 1e-10 ){

    J1 = which.min(abs(b))  # establish the zero direction
    J2 = (J1 %% 3) + 1;
    J3 = ((J1+1) %% 3)+1;
    A1 = a[J1]
    A2 = a[J2]
    A3 = a[J3]
    B1 = b[J1]
    B2 = b[J2]
    B3 = b[J3]

    # Calculate polynomial co-efficients
    coeff = c( 2*(B2^2+B3^2),
      B2^2*(-2*A2-5*A3-A1)/10+B3^2*(-2*A3-5*A2-A1)/10,
      2*B2^2*(2*A3^2+2*A2*A3+A3*A1)/100+2*B3^2*(2*A2^2+2*A2*A3+A2*A1)/100,
      B2^2*(-2*A2*A3^2-A3^3-A3^2*A1)/1000+B3^2*(-2*A3*A2^2-A2^3-A2^2*A1)/1000)

    # Solve polynomial
    mus = polyroot(coeff[4:1]);

    # Find horizon
    horizon = matrix(rep(0,3))
    for(k in 1:length(mus)){
      
      if ( abs(Im(mus[k] )) > 1e-6){next }    # Skip solutions with imaginary component
      
      nr_test = c(0,0,0)
      nr_test[J2] = B2/(10*Re(mus[k])-A2)/2   # Calculate horizon
      nr_test[J3] = B3/(10*Re(mus[k])-A3)/2

      if(max(is.na(nr_test) | is.infinite(nr_test))){next}   # Skip if there are NaNs or infs

      # Only keep if bigger than previous horizon and is inside the Bloch ball
      if ( norm(nr_test, type="2") > norm(horizon, type="2") & norm(nr_test, type="2") <= 1){
        horizon = nr_test
      }

    }
    return(horizon)            
  }
  
  # For b with all non-zero elements, we must solve a quintic equation
  coeff_exp = matrix(0,6,3);
  for(j in 1:3){
    j1 = (j %% 3) + 1;
    j2 = ((j+1) %% 3)+1;
    A1 = a[j]; B1 = b[j];
    A2 = a[j1]; B2 = b[j1];
    A3 = a[j2]; B3 = b[j2];
  
    # Calculate polynomial  co-efficients
    coeff_exp[1,j] = 2;   
    coeff_exp[2,j] = -4*(A2+A3);   
    coeff_exp[3,j] = 2*(A2^2+A3^2+4*A2*A3);   
    coeff_exp[4,j] = -4*(A2*A3^2+A3*A2^2);   
    coeff_exp[5,j] = 2*A2^2*A3^2;   
    coeff_exp[2:6,j] = coeff_exp[2:6,j] - coeff_exp[1:5,j]*(A2/2+A3/2+A1);
    coeff_exp[,j] = coeff_exp[,j]*B1^2;
  }
  coeff = rowSums(coeff_exp);

  # Solve polynomial  
  mus = polyroot(coeff[6:1]);

  # Calculate horizon
  horizon = matrix(rep(0,3));
  for(k in 1:length(mus)){
    if ( abs(Im(mus[k] )) > 1e-3){next } # Skip complex solutions
    nr_test = b/(Re(mus[k])-a)/2         # Calculate horizon
    if(max(is.na(nr_test) | is.infinite(nr_test))){next} # Skip Nans and infs
  
    # Keep only if new horizon is bigger than old, and is inside the Bloch ball
    if ( norm(nr_test, type="2") > norm(horizon, type="2") & norm(nr_test, type="2") <= 1){
      horizon = nr_test}
  }
  return(horizon)
}