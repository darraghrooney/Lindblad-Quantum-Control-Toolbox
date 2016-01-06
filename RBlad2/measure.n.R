# List of functions in this file:
#
#   measure.n()
#   xyz.measure.n()
#   rand.measure.n.unn()
#   grid.measure.n.unn()
#   rand.n()
#   measurements.generate()
#   create.training.mm()
#   create.crossval.mm()
#   create.test.mm()

# Functions required for many in this file:
source("system.functions.R")

# This function simulates a measurement of n given its actual value
# The number of trials is given as an input variable. The measurements are made along
# the x,y,z axes.

measure.n <- function(no.meas, actual.n, direction){
  direction = direction / sqrt(sum(direction^2))
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

rand.measure.n.unn <- function(n, mm){
  n.rand = c(0,0,0)
  for (k in 1:mm){
    m = rand.n()
    m = m/mag(m)
    n.rand = n.rand + measure.n(1, n, m)*m
  }
  return(n.rand)  
}

grid.measure.n.unn <- function(n, mm){
  sqnm = floor(sqrt(mm))
  n.grid = c(0,0,0)
  for (k in 1:mm){
    cosphi = 2*( (k-1) / (mm-1) ) -1
    theta = floor( (k-1) %% sqnm)*2*pi/sqnm
    m = c(sqrt(1-cosphi^2)*c(cos(theta),sin(theta)),cosphi)
    n.grid = n.grid + measure.n(1,n,m)*m
  }
  return(n.grid)
}


# Returns a random vector in the Bloch ball. Longitude and latitude are uniform,
# and the cube root of the magnitude is uniform. This gives a uniform distribution
# on the ball.

rand.n <- function(){
  nmag = runif(1)^(1/3)
  cosphi = 2*runif(1)-1
  theta = runif(1)*2*pi
  return(nmag*c(sqrt(1-cosphi^2)*c(cos(theta),sin(theta)), cosphi))
}


# This is a function for compiling data on measuring states in the Bloch ball.
# There are three methods: (1) measuring along the x-y-z axes (which are determined a
# priori), (2) measuring along random axes and (3) measuring along axes determined by
# a spiral-grid method. The spiral-grid method consists of linearly incrementing the
# z-coordinate and the theta-coordinate, the ratio being ~ the square-root of the 
# data size. 

# The data that is compiled: (1) the actual Bloch vector and its magnitude, (2) its estimation
# and magnitude using the xyz method, as well as the magnitude/angle/total errors 
# (3) the angular estimation using the random method, together with the compiled
# magnitude, ratio of estimated to actual magnitude, and angular error, and
# (4) everything in (3), but using the spiral-grid method.

# The actual Bloch vectors are random, in the sense that the distribution is uniform
# over the ball. The first input is the number of Bloch vectors, the second a vector of 
# the number of measurements for each.

# meas is a vector of measurement numbers. npm specifies how many random Bloch vectors are to
# be measured that number of times. The seed can also be set.

measurements.generate <- function(meas, npm  = 1, seed = 53){

  set.seed(seed)
  
  # Compute size of data frame and intialize

  meas.no = length(meas)
  nr = npm*meas.no
  results = data.frame( meas = rep(meas, each = npm),
        n1 = double(nr), n2 = double(nr), n3 = double(nr), nmag = double(nr),
        xyz.n1 = double(nr), xyz.n2 = double(nr), xyz.n3 = double(nr), xyz.nm = double(nr),
        xyz.mag.err = double(nr), xyz.ang.err = double(nr), xyz.tot.err = double(nr),
        rand.n1 = double(nr), rand.n2 = double(nr), rand.n3 = double(nr), 
        rand.nm = double(nr), rand.mr = double(nr), rand.ang.err = double(nr), 
        grid.n1 = double(nr), grid.n2 = double(nr), grid.n3 = double(nr),
        grid.nm = double(nr), grid.mr = double(nr), grid.ang.err = double(nr))
  
  # Loop over measurement numbers
  
  for (l in 1:meas.no){
  
    print(paste("Measurement", l, "out of", meas.no, sep=" "))
    mm = meas[l]

    # Loop over vectors
    for (j in 1:npm){

      print(paste("... vector ", j, ".", l, sep=""))
      
      # Generate actual Bloch vector
      n = rand.n()
      n.mag = mag(n)

      # Measure using the x-y-z method
      n.xyz = xyz.measure.n(n, mm)
      n.xyz.mag = mag(n.xyz)

      # Calculate errors in the x-y-z method
      n.xyz.me = n.xyz.mag - n.mag
      n.xyz.ae = 1-(n %*% n.xyz)/ n.mag / n.xyz.mag
      n.xyz.te = mag(n-n.xyz)
      
      # Measure using the (unnormalized) random-axis method
      n.rand = rand.measure.n.unn(n, mm)
      n.rand.mag = mag(n.rand)
      
      # Compute angular err
      n.rand.ae = 1 - (n %*% n.rand)/ n.mag / n.rand.mag
    
      # Measure using spiral-grid method
      sqnm = floor(sqrt(mm))
      n.grid = c(0,0,0)
      for (k in 1:mm){
        cosphi = 2*( (k-1) / (mm-1) ) -1
        theta = floor( (k-1) %% sqnm)*2*pi/sqnm
        m = c(sqrt(1-cosphi^2)*c(cos(theta),sin(theta)),cosphi)
        n.grid = n.grid + measure.n(1,n,m)*m
      }
      n.grid.mag = mag(n.grid)
      
      # Compute angular error
      n.grid.ae = 1 - (n %*% n.grid)/ n.mag / n.grid.mag

      # Record results      
      results[(l-1)*npm + j, -1] = c(n, n.mag, n.xyz, n.xyz.mag, n.xyz.me, n.xyz.ae, n.xyz.te, 
                n.rand, n.rand.mag, n.rand.mag/n.mag, n.rand.ae, 
                n.grid, n.grid.mag, n.grid.mag/n.mag, n.grid.ae)
    }
  }    
  return(results)
}

# The following three functions construct data sets for the measurement portion of the project. The
# training set has 100 different measurement numbers between 30 and 3000, and there are 700 Bloch vectors
# for each. The cross-val set has 30000 randomized meas nos between 30 and 6000, and the test set has 
# 100000. The meas nos are not uniformly distributed: smaller ones are weighted (the sqaures are uniformly
# distributed).

create.training.mm <- function(){

	meas.lims = c(30,3000)
  	meas.no = 100
  	mms = round((min(sqrt(meas.lims))+diff(range(sqrt(meas.lims)))*(0:(meas.no-1))/(meas.no-1))^2)
  	return(measurements.generate(mms, 700, 53))
}

create.crossval.mm <- function(){
  	n = 30000
  	meas.lims = c(30,6000)
  	mms = round((min(sqrt(meas.lims))+diff(range(sqrt(meas.lims)))*runif(n))^2)
  	return(measurements.generate(mms, seed = 151))  
}

create.test.mm <- function(){
  	n = 100000
  	meas.lims = c(30,6000)
  	mms = round((min(sqrt(meas.lims))+diff(range(sqrt(meas.lims)))*runif(n))^2)
  	return(measurements.generate(mms, seed = 23432))  
}