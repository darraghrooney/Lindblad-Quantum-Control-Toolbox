# This file contains the following functions:
#
# ellipsoid.plot.R
# goldilocks.example.R

# Load pre-requisite code
source("algorithm.R")
source("sys.solve.R")
library(scatterplot3d)

# This function plots an example ellipsoid attractor (in blue) encased in the 
# Bloch sphere (in grey). The "naked" switch turns off the box, grid and axes.

ellipsoid.plot <- function(a, b, naked = FALSE){ # We assume no rotation
  
  # color vectors
  
  incol = rgb(0,0,1)
  outcol = rgb(0,.85,0)
  
  # Construct Bloch sphere points
  bb.grid = 100
  bb.th = (1:bb.grid)/bb.grid*2*pi
  bb.cp = 2*(0:bb.grid)/bb.grid - 1
  bb.sp = sqrt(1-bb.cp^2)
  bb = data.frame(nx = kronecker(cos(bb.th), bb.sp), 
                  ny = kronecker(sin(bb.th), bb.sp),
                  nz = kronecker(rep(1,bb.grid), bb.cp))
  
  # Plot the back of the Bloch sphere.
  if (naked){
    s3d <- scatterplot3d(bb[1:((bb.grid+1)*bb.grid*.6),], cex.symbols = .4, pch = 18, 
                         color = outcol, angle = 10, box=FALSE, grid = FALSE, axis=FALSE)
  }
  else{
    s3d <- scatterplot3d(bb[1:((bb.grid+1)*bb.grid*.6),], cex.symbols = .4, pch = 18, 
                         color = outcol, angle = 10)
  }
  
  # Construct ellipsoid
  a.comp = sum(a) - a
  centre = b/2/a.comp
  axes = sqrt(sum(a.comp*centre*centre)/a.comp)
  am = data.frame(nx = centre[1] + axes[1]*kronecker(cos(bb.th), bb.sp), 
                  ny = centre[2] + axes[2]*kronecker(sin(bb.th), bb.sp),
                  nz = centre[3] + axes[3]*kronecker(rep(1,bb.grid), bb.cp))

  # Plot ellipsoid
  s3d$points(am,cex =.3, pch=18, col=incol)
  
  # Plot front of Bloch sphere
  s3d$points(bb[((bb.grid+1)*bb.grid*.6):((bb.grid+1)*bb.grid),], 
             cex = .4, pch = 18, col = outcol)  
}

# This function is used to illustrate the problem of scale: parameter estimation
# works if our assumed scale is just right, but horrendous if too big or small
# The input H indicates our presumed scale. The "plt" switch plots the collected
# measurements

goldilocks.example <- function(H, system, plt = TRUE){

  # Get measurements
  mm = collect.mm(system, mph = 1e4, H=H, naive = TRUE)
  
  # Plot if desired (if not, return the measurements)
  if (plt){
    x = calc.x(as.matrix(mm[,1:3]), as.matrix(mm[,7:9]))
    names(mm)[7:9]=c("nx","ny", "nz")

    # Plot
    scatterplot3d(mm[,7:9], xlim = c(-1,1), ylim=c(-1,1), zlim=c(-1,1),cex.symbols=.4,
                  pch=18)
    
    # Print table comparing actual and estimated parameters
    comparison = data.frame(cbind(c(system$b, system$a,0,0,0),x,abs(x-c(system$b, system$a,0,0,0))))
    row.names(comparison) = c("b1","b2","b3","a11","a22","a33","a23","a31","a12")
    names(comparison) = c("Actual value", "Estimated value", "Discrepancy")
    print(comparison)
  }
  else{
    return(mm)
  }
}