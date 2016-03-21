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
