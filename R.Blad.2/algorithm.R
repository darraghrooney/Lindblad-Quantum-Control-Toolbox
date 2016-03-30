# This file contains the following functions:
#
# collect.mm()
# system.determine()
# estimate.system()

# Pre-reqs
source("system.functions.R")
source("measure.n.R")
source("sys.solve.R")

# This function takes measurements of a system. The hyper-parameters mph, p1 and p2 
# must be specified. The "naive" switch turns off the scale tuning
collect.mm <- function(system, mph, p1, p2, naive = FALSE){
  
  # Parameters for assigning measurements
  ham.no = ceiling(MM.BUDGET/mph)
  leftover = MM.BUDGET - mph*(ham.no-1)
  
  # Initialize scale to least possible a_1 scale
  est.scale = MIN.SCALE
  
  a = system$a
  b = system$b
  rot = rot.from.u(system$u)
  
  # Parameters to construct Hamiltonians
  thetas = runif(ham.no-1)*2*pi
  cosphis = 2*runif(ham.no-1)-1
  sinphis = sqrt(1-cosphis^2)
  h.mags = runif(ham.no-1)
  
  measurements = data.frame(hx = c(0,h.mags*sinphis*cos(thetas)), 
                            hy = c(0,h.mags*sinphis*sin(thetas)), hz = c(0,h.mags*cosphis), 
                            nx.actual = numeric(ham.no), ny.actual = numeric(ham.no), 
                            nz.actual = numeric(ham.no), nx.meas = numeric(ham.no), 
                            ny.meas = numeric(ham.no), nz.meas = numeric(ham.no))
  
  # Sweep over hamiltonians
  for (j in 1:ham.no){
    
    # Construct hamiltonian
    ham = est.scale*measurements[j, 1:3]
    measurements[j, 1:3] = ham
    
    # Fetch actual Bloch vector
    n.actual = as.vector(rot %*% n.from.h(t(rot) %*% 
                                            t(as.matrix(ham)),a,b))
    measurements[j, 4:6] = n.actual
    
    # Measure Bloch vector
    if (j < ham.no){
      n.meas = xyz.measure.n(n.actual, mph)
    }
    else{
      n.meas = xyz.measure.n(n.actual,leftover)
    }
    measurements[j, 7:9] = n.meas
    
    # If not naive, adjust estimated scale according to mean and sd
    # of as-yet-measured Bloch vectors
    if(!naive){
      est.mean = mean(mag(measurements[1:j,7:9]))
      est.sd = 0
      if (j>1){
        est.sd = sd(mag(measurements[1:j,7:9]))
      }
      if (sum(is.na(c(est.mean,est.sd,p2))))
        {print(c(est.mean,est.sd,p2))}
    
      if (est.sd < p2*est.mean){
        est.scale = est.scale*p1^(1-est.sd/p2/est.mean)
      }
    } 
  }
  return(measurements)
}

# This function determines the parameters of the system. In god-mode, uses
# actual Bloch vector. In ungod-mode, measures the Bloch vector and estimates.

system.determine <- function(system, mph, p1, p2, 
                             god = FALSE, naive = FALSE){

  # Fetches the actual and measured Bloch vectors
  mms = collect.mm(system, mph, p1, p2, naive)
  ham.no = dim(mms)[1]
  
  # Get actual parameters
  if (god){
    sys.est = calc.x(mms[,1:3],mms[,4:6])
  }
  # Estimate parameters
  else{
    sys.est = calc.x(mms[,c(1:3)], mms[,c(7:9)] )
  }      
  return(sys.est)
}

