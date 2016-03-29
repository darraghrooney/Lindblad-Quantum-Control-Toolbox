# List of functions in this file:
#
#   measure.n()
#   xyz.measure.n()

# Functions required for many in this file:
source("system.functions.R")

# This function simulates a measurement of n given its actual value
# The number of trials is given as an input variable. The measurements are made along
# the x,y,z axes.

measure.n <- function(no.meas, actual.n, direction){
  direction = direction / sqrt(sum(direction^2))
  if(class(actual.n) == "data.frame"){
    actual.n = as.matrix(actual.n)
  }
  prob = (actual.n %*% direction + 1)/2        # find probabilities
  return( 2*rbinom(1, no.meas, prob)/no.meas - 1) # generate measurement
}

# The following three functions are various measurement methods for Bloch vectors.
# See comments for measurements.generate(). The latter two yield only the angular
# component, as they require a posteriori multiplicative constant.

xyz.measure.n <- function(n, meas.no){
  return(c( measure.n(floor(meas.no/3), n, c(1,0,0)), 
            measure.n(floor(meas.no/3), n, c(0,1,0)),
            measure.n(meas.no - 2*floor(meas.no/3), n, c(0,0,1))))
}

collect.mm <- function(system, mph, p1 = 1, p2 = 0, H = 1e4, naive = TRUE){
  
  ham.no = ceiling(MM.BUDGET/mph)
  leftover = MM.BUDGET - mph*(ham.no-1)
  est.scale = H
  
  a = system$a
  b = system$b
  rot = rot.from.u(system$u)
  
  thetas = runif(ham.no-1)*2*pi
  cosphis = 2*runif(ham.no-1)-1
  sinphis = sqrt(1-cosphis^2)
  h.mags = runif(ham.no-1)
  
  measurements = data.frame(hx = c(0,h.mags*sinphis*cos(thetas)), 
            hy = c(0,h.mags*sinphis*sin(thetas)), hz = c(0,h.mags*cosphis), 
            nx.actual = numeric(ham.no), ny.actual = numeric(ham.no), 
            nz.actual = numeric(ham.no), nx.meas = numeric(ham.no), 
            ny.meas = numeric(ham.no), nz.meas = numeric(ham.no))
  
  mean.sum = 0
  var.sum = 0
  for (j in 1:ham.no){
    ham = est.scale*measurements[j, 1:3]
    measurements[j, 1:3] = ham
    n.actual = as.vector(rot %*% n.from.h(t(rot) %*% 
                                            t(as.matrix(ham)),a,b))
    if (j < ham.no){
      n.meas = xyz.measure.n(n.actual, mph)
    }
    else{
      n.meas = xyz.measure.n(n.actual,leftover)
    }
    measurements[j, 4:6] = n.actual
    measurements[j, 7:9] = n.meas
    
    if(!naive){
      mean.sum = mean.sum + mag(n.meas)
      est.mean = mean.sum/j
      var.sum = var.sum + (mag(n.meas)-est.mean)^2
      est.sd = 0
      if (j>1){
        est.sd = sqrt(var.sum/(j-1))
      }
      if (est.sd < p2*est.mean){
        est.scale = est.scale*p1^(1-est.sd/p2/est.mean)
      }
    } 
  }
  return(measurements)
}
