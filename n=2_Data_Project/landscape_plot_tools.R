# These are functions to produce various plots for the project

# Plots the discretization of a-space. Can plot points or lines.
plot_as <- function(rot_ang=70, pt="l"){

  # Generate data
  x=rep(100,26*27/2)
  y=c(0)
  z=c(0)
  for (j1 in 1:25){
    for (j2 in 0:j1){
      y=c(y,4*j1)
      z=c(z,4*j2)
    }
  }
  coord=data.frame(x, y, z)
  
  # Create plot object
  library(scatterplot3d)
  sp3 = scatterplot3d(coord$x,coord$y,coord$z, type="n", box=FALSE,axis=FALSE,
                      grid=FALSE, angle=rot_ang, scale.y=.8)
  
  # The parameter space is divided in 4 "stripes", each with a different color
  mycols = c("magenta", "blue","green", "red")
  for (j in 0:25){
    line_to_plot = coord[coord$z==4*j,]
    sp3$points3d(line_to_plot$x,line_to_plot$y,line_to_plot$z, type=pt, 
                 col=mycols[1+4*j/26], pch=19, cex=.5)
  }

  # The boundary points are given a solid black line
  line_to_plot = coord[coord$y==coord$z,]
  sp3$points3d(line_to_plot$x,line_to_plot$y,line_to_plot$z, type="l", 
               col="black", lwd=2)
  line_to_plot = coord[coord$y==100,]
  sp3$points3d(line_to_plot$x,line_to_plot$y,line_to_plot$z, type="l", 
               col="black", lwd=2)
  line_to_plot = coord[coord$z==0,]
  sp3$points3d(line_to_plot$x,line_to_plot$y,line_to_plot$z, type="l", 
               col="black", lwd=2)
  }

# Plots the discretization of angle-space
plot_angles <- function(rot_ang=70, pt="l"){

  # Generate data
  cp = rep(0:30, 31)
  th = rep(0:30, 1, each=31)
  
  z=cp/30
  x=1-sqrt(1-(cp/30)^2)*cos(th/30*pi/2)
  y=1-sqrt(1-(cp/30)^2)*sin(th/30*pi/2)

  coord=data.frame(cp, th, x, y, z)
    
  # Create plot object
  library(scatterplot3d)
  sp3 = scatterplot3d(coord$x,coord$y,coord$z, type="n", box=FALSE, axis=FALSE,
        grid=FALSE, angle=rot_ang, scale.y=0.9)
  
  # Angle-space is an eighth of a sphere. This space is also striped.
  mycols = c("magenta", "blue","green", "red")
  for (j in 0:30){
    line_to_plot = coord[coord$cp==j,]
    sp3$points3d(line_to_plot$x,line_to_plot$y,line_to_plot$z, type=pt, 
                 col=mycols[1+4*j/31], pch=19, cex=.5)
  }

  # Solid boundary line
  line_to_plot = coord[coord$cp==0,]
  sp3$points3d(line_to_plot$x,line_to_plot$y,line_to_plot$z, type="l", 
               col="black", lwd=2)
  line_to_plot = coord[coord$th==0,]
  sp3$points3d(line_to_plot$x,line_to_plot$y,line_to_plot$z, type="l", 
               col="black", lwd=2)
  line_to_plot = coord[coord$th==30,]
  sp3$points3d(line_to_plot$x,line_to_plot$y,line_to_plot$z, type="l", 
               col="black", lwd=2)

  # The sphere axes are also drawn in with dashed lines for perspective
  line_to_plot = data.frame(x=rep(1,6),y=rep(1,6),z=(0:5)/5)
  sp3$points3d(line_to_plot$x,line_to_plot$y,line_to_plot$z, 
               type="l", col="black", lty=2)
  line_to_plot = data.frame(x=rep(1,6),z=rep(0,6),y=(0:5)/5)
  sp3$points3d(line_to_plot$x,line_to_plot$y,line_to_plot$z, 
               type="l", col="black", lty=2)
  line_to_plot = data.frame(z=rep(0,6),y=rep(1,6),x=(0:5)/5)
  sp3$points3d(line_to_plot$x,line_to_plot$y,line_to_plot$z, 
               type="l", col="black", lty=2)
}

# This functions plots three variables (out of a possible ten) while
# fixing either a2_ind and a3_ind, or cosphi_ind and theta_ind, and
# sweeping over the complementary pair. vars_to_fix should be either 
# c(1,2) or c(3,4), and vars_to_plot should be a triple from 1:10
# (that is disjoint from vars_to_fix). The plots are striped, so 
# Var_to_stripe should be in 1:4 but not in vars_to_fix

plot_fix_var <- function(vars_to_fix, fixed_values, vars_to_plot,
    var_to_stripe, pt = "p", rot_ang=150, need_to_read=FALSE){
  
  # Read in landscape data if necessary
  if (need_to_read){
    library(R.matlab)
    rawls = readMat("landscape.MAT")
    landscape = data.frame(rawls$landscape[,,1])
  }
  
  library(scatterplot3d)
  
  # Select out the data for fixed values
  to_plot = landscape[landscape[,vars_to_fix[1]]==fixed_values[1] &
        landscape[,vars_to_fix[2]]==fixed_values[2],]
  axlab = names(landscape)[vars_to_plot]
  
  # Create plot object
  sp3 = scatterplot3d(to_plot[,vars_to_plot[1]],to_plot[,vars_to_plot[2]],
        to_plot[,vars_to_plot[3]], type="n", grid=FALSE, angle=rot_ang, 
        xlab=axlab[1], ylab=axlab[2], zlab=axlab[3], box=FALSE)
  
  mycols = c("magenta", "blue","green", "red")
  jmax = max(to_plot[,var_to_stripe])
  
  # Sweep over data, and stripe it
  for (j in 0:jmax){
    line_to_plot = to_plot[to_plot[,var_to_stripe]==j,]
    sp3$points3d(line_to_plot[,vars_to_plot[1]],line_to_plot[,vars_to_plot[2]],
          line_to_plot[,vars_to_plot[3]], type=pt, cex=.5, pch=16,
          col=mycols[1 + 4*j/(jmax+1)] )
  }
  
  # Plot boundary lines
    if(length(setdiff(c(1,2),vars_to_fix))==0){
      line_to_plot = to_plot[to_plot$theta_ind==0,]
      sp3$points3d(line_to_plot[,vars_to_plot[1]],
          line_to_plot[,vars_to_plot[2]],line_to_plot[,vars_to_plot[3]], 
                   type="l", col="black", lwd=2)
      line_to_plot = to_plot[to_plot$theta_ind==30,]
      sp3$points3d(line_to_plot[,vars_to_plot[1]],
                   line_to_plot[,vars_to_plot[2]],line_to_plot[,vars_to_plot[3]], 
                   type="l", col="black", lwd=2)
      line_to_plot = to_plot[to_plot$cosphi_ind==0,]
      sp3$points3d(line_to_plot[,vars_to_plot[1]],
                   line_to_plot[,vars_to_plot[2]],line_to_plot[,vars_to_plot[3]], 
                   type="l", col="black", lwd=2)
    }
    else if(length(setdiff(c(3,4),vars_to_fix))==0){
        line_to_plot = to_plot[to_plot$a3_ind==0,]
        sp3$points3d(line_to_plot[,vars_to_plot[1]],
                     line_to_plot[,vars_to_plot[2]],line_to_plot[,vars_to_plot[3]], 
                     type="l", col="black", lwd=2)
        line_to_plot = to_plot[to_plot$a2_ind==25,]
        sp3$points3d(line_to_plot[,vars_to_plot[1]],
                     line_to_plot[,vars_to_plot[2]],line_to_plot[,vars_to_plot[3]], 
                     type="l", col="black", lwd=2)
        line_to_plot = to_plot[to_plot$a2_ind==to_plot$a3_ind,]
        sp3$points3d(line_to_plot[,vars_to_plot[1]],
                     line_to_plot[,vars_to_plot[2]],line_to_plot[,vars_to_plot[3]], 
                     type="l", col="black", lwd=2)
      }
}
