# This file contains the following files:
#
# make.training.set()
# make.test.set()

# Pre-reqs
source("system.functions.R")
source("measure.n.R")
source("algorithm.R")
library(R.matlab)

# This function makes the training set. There are 7000 randomized systems. The 
# hyper-parameters p1, p2 and M are chosen randomly, as well as the scale.

make.training.set <- function(seed, dest.file, po.no){
  set.seed(281)
  po.no = 7000
  train.set = data.frame(matrix(0,po.no,19))
  names(train.set) = c("scale","meas.no", "p1","p2","a1.act","a2.act",
        "a3.act","b1.act","b2.act","b3.act","a11.meas",
        "a22.meas","a33.meas","a12.meas",
        "a23.meas","a31.meas","b1.meas","b2.meas","b3.meas")

  # Sweep over systems
  for (j in 1:po.no){
    print(paste("Estimating system", j))
    
    # Set scale and generate system
    scale = MIN.SCALE*(MAX.SCALE/MIN.SCALE)^runif(1)
    system = rand.A(scale, FALSE)

    # Generate p1, p2 and meas (which is MM.BUDGET/M)
    p1 = 1+2*runif(1)
    p2 = runif(1)/2
    meas = 1e4+round(9e4*runif(1))
    
    # Estimate parameters and add to data set
    estimate = system.determine(system, meas, p1, p2) 
    train.set[j,] = c(scale, meas, p1, p2, system$a, system$b, estimate[4:9],estimate[1:3])
    
  }
  writeMat("train.set.MAT", training = train.set)
}

# Generates test set. The difference is that the hyper-parameters are now set.

make.test.set <- function(){
  
  # Fetch p1, p2
  p.opt = readMat("p.opt.MAT")$p.opt[]
  
  set.seed(280)
  test.size = 3000

  # Set hyper-parameters
  mph = 5e4
  p1 = p.opt[1]
  p2 = p.opt[2]
  test.set = data.frame(matrix(0, test.size, 15))
  names(test.set) = c("a1.act","a2.act", "a3.act","b1.act","b2.act","b3.act","a11.meas",
                      "a22.meas","a33.meas","a12.meas",
                      "a23.meas","a31.meas","b1.meas","b2.meas","b3.meas")
  
  for (j in 1:test.size){
    print(paste("Estimating system", j))
    scale = MIN.SCALE*(MAX.SCALE/MIN.SCALE)^runif(1)
    system = rand.A(scale, FALSE)
    estimate = system.determine(system, mph, p1, p2) 
    test.set[j,] = c(system$a, system$b, estimate[4:9],estimate[1:3])
  }
  
  writeMat("test.set.MAT", test = test.set)
}