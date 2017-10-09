# This file contains the following files:
#
# p.optimize()
# best.pt()
# circumradius()

# Pre-reqs:
library(deldir)
library(R.matlab)
source("assessment.R")
train.set = data.frame(readMat("assess.train.MAT")$assess[,,1])

# This function optimizes p1 and p2 by computing the Delaunay triangulation
# and then using the point with largest circumradius

p.optimize <- function(){

  train.size = dim(train.set)[1]

  # Any estimation that has a metric worse than 1-exp(-2) ~ 0.86 is considered
  # "bad" 
  tol = -2

  # Find nans and bad points
  nan.pts = is.nan(train.set$trace.metric) & is.nan(train.set$diag.metric)
  bad.pts1 = log(1-train.set$trace.metric) > tol
  bad.pts1[is.na(bad.pts1)] = TRUE
  bad.pts2 = log(1-train.set$diag.metric) > tol
  bad.pts3 = log(1-train.set$off.diag.metric) > tol
  bad.pts4 = log(1-train.set$b.metric) > tol
  bad.set = train.set[nan.pts | bad.pts1 | bad.pts2 | bad.pts3 | bad.pts4,]
  bad.size = dim(bad.set)[1]

  # Compute triangulation and plot the Dirichlet tessellation  
  dd = deldir(bad.set[,c("p1", "p2")])
  plot.tile.list(tile.list(dd),border = 2,xlim=c(1,3),ylim=c(0,.5))
  
  # Find the optimal p1, p2
  p.opt = best.pt(dd)
  
  # Return both the optimization and the set of bad points
  return(list(bad.set = bad.set, p.opt = p.opt))
}

# This function takes a Delaunay triangulation and finds the optimal point 
# (essentially the triangle with the largest circumradius)

best.pt <- function(dd){
  
  # Get the triangles
  deltri = triang.list(dd)
  best.pt = c(0,0,0)
  
  # Sweep over triangles to get the optimal points
  for (j in 1:length(deltri)){
    new.pt = circumradius(deltri[[j]])
    
    # Ignore points outside the allowed parameter space
    if (new.pt$cr > best.pt[3] & new.pt$cc[1] > 1  & new.pt$cc[1] < 3
        & new.pt$cc[2] > 0 & new.pt$cc[2] < 0.5){
      best.pt = c(new.pt$cc, new.pt$cr)
    }
  }
  return(best.pt)
}

# This function computes the circumradius and circumcenter of a triangle
circumradius <- function(triangle){
  xs = triangle$x
  ys = triangle$y
  msqs = xs^2 + ys^2
  a = det(matrix(c(xs,ys,1,1,1),3,3))
  b = det(matrix(c(xs,ys,msqs),3,3))
  Sx = .5*det(matrix(c(msqs,ys,1,1,1),3,3))
  Sy = .5*det(matrix(c(xs,msqs,1,1,1),3,3))

  return(list(cc = c(Sx/a,Sy/a), cr = sqrt(b/a+(Sx^2+Sy^2)/a^2)))
}
