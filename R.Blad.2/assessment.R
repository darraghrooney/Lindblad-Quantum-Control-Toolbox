# This file contains the following functions
#
# assess()
# assess.test()

# Pre-reqs
library(R.matlab)
source("system.functions.R")

# This function takes the training data and calculates four different metrics for each
# parameter estimation: one for scale (aka trace), one for the diagonals of A (independent of 
# trace), one for the off-diagonals of A, and one for b

assess <- function(){

  # Get training data and initialize assessment
  df = data.frame(readMat("train.set.MAT")$training[,,1])
  assessment = data.frame(matrix(0,dim(df)[1],8))
  names(assessment) = c("scale", "meas.no", "p1", "p2", 
    "trace.metric", "off.diag.metric","diag.metric", "b.metric")

  # Loop over estimationes
  for (j in 1:dim(df)[1]){
    real.params = df[j,5:10]
    est.params = df[j,11:19]
    real.scale = sum(real.params[1:3])/3
    
    # Calculate metrics
    trace.metric = 1 - abs(log(sum(est.params[1:3])/3/real.scale))
    off.diag.metric = 1 - mag(est.params[4:6])/real.scale
    diag.metric = 1 - sqrt( ((real.params[1]-real.params[2]) - 
                               (est.params[1]-est.params[2])  )^2/2 + 
                              ((real.params[2]-real.params[3]) - 
                                 (est.params[2]-est.params[3])   )^2/2 ) / real.scale
    b.metric = 1- mag(real.params[4:6] - est.params[7:9])/real.scale
    assessment[j,] = c(df[j,1:4],trace.metric,off.diag.metric,diag.metric,
                       b.metric)
  }
  
  return(assessment)
}

# The same function as above, but for the test set. The returned df omits the 
# hyper-parameters scale, p1, p2 and meas.no
assess.test <- function(){
  
  df = data.frame(readMat("test.set.MAT")$test[,,1])
  assessment = data.frame(matrix(0,dim(df)[1],4))
  names(assessment) = c("trace.metric", "off.diag.metric","diag.metric", "b.metric")
  for (j in 1:dim(df)[1]){
    print(j)
    real.params = df[j,1:6]
    est.params = df[j,7:15]
    real.scale = sum(real.params[1:3])/3
    trace.metric = 1 - abs(log(sum(est.params[1:3])/3/real.scale))
    off.diag.metric = 1 - mag(est.params[4:6])/real.scale
    diag.metric = 1 - sqrt( ((real.params[1]-real.params[2]) - 
          (est.params[1]-est.params[2])  )^2 + 
          ((real.params[2]-real.params[3]) - 
          (est.params[2]-est.params[3])   )^2 ) / real.scale
    b.metric = 1- mag(real.params[4:6] - est.params[7:9])/real.scale
    assessment[j,] = c(trace.metric,off.diag.metric,diag.metric,
                       b.metric)
  }
  
  return(assessment)
}
