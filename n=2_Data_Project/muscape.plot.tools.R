# These are functions to produce various plots for the project

# Plots the discretization of a-space.
plot.as <- function(ms.data){

  # Create data frame
  a.data = data.frame(unique(ms.data[,1:2]))
  nr = dim(a.data)[1]
  a.data = data.frame(a.data, color = character(nr))

  # Assign colors by striping according to the a3 variable
  mycols = c("magenta", "blue","green", "red", "darkorange", "violet")
  a.data$color = mycols[findInterval(a.data$a3, c(0,20,40,60,80,100))]

  # Plot
  plot(c(0,100,100,0), c(0,0,100,0), xlab="a2", ylab="a3", type="l")  # Boundary
  points(a.data$a2,a.data$a3, col = a.data$color, pch = 19, cex = .5) # Data
  par(mar=c(5,4,4,2))   # I needed to adjust left margin to see the y-axis label
}

# Plots the discretization of angle-space
plot.angles <- function(ms.data, rot.ang = 70){
  
  # Create data frame
  angle.data = data.frame(unique(ms.data[,3:4]))
  nr = dim(angle.data)[1]
  sp = sqrt(1-angle.data$cosphi^2)
  angle.data = data.frame(angle.data, bx = sp*cos(angle.data$theta), 
                  by = sp*sin(angle.data$theta), color = character(nr))
  
  # Assign colors according to cosphi variable
  mycols = c("magenta", "blue","green", "red", "darkorange", "violet")
  angle.data$color = mycols[findInterval(angle.data$cosphi, 0:5/4.999)]
  
  # Plot
  library(scatterplot3d)
  ad = 0:99/200*pi
  sp3 = scatterplot3d(angle.data$bx,angle.data$by, angle.data$cosphi,  
            grid=FALSE, box = FALSE, type="n", angle = rot.ang, xlab = "b1",
            ylab="b2", zlab="b3")
  # Boundary
  sp3$points(c(cos(ad),rep(0,99),sin(ad)), c(sin(ad),cos(ad),rep(0,99)), 
                      c(rep(0,99),sin(ad),cos(ad)), type="l")
  # Data
  sp3$points(angle.data$bx,angle.data$by, angle.data$cosphi, col = angle.data$color, 
         pch = 19, cex = .5)
}

# This functions plots three variables (out of a possible twelve) while
# fixing either a2 and a3, or cosphi and theta, and
# sweeping over the complementary pair. vars.to.fix should be either 
# c(1,2) or c(3,4), and vars.to.plot should be a triple from 1:12
# (that is disjoint from vars.to.fix). The plots are striped, so 
# var.to.stripe should be in 1:4 but not in vars.to.fix

plot.fix.2var <- function(ms.data, vars.to.fix, fixed.values, vars.to.plot,
    var.to.stripe, rot.ang=150, a.reflect=TRUE){
  
  library(scatterplot3d)
  
  # Select out the data for fixed values
  eps = 1e-8
  to.plot = ms.data[ findInterval( ms.data[,vars.to.fix[1]], 
        fixed.values[1] + c(-eps,eps)) == 1  &
        findInterval( ms.data[,vars.to.fix[2]], fixed.values[2]+c(-eps,eps))==1, 
        c(setdiff(1:4,vars.to.fix), vars.to.plot,var.to.stripe) ]
 
  # Assign colors
  colormap = character(dim(to.plot)[1])
  mycols = c("magenta", "blue","green", "red", "darkorange", "violet")
  if (var.to.stripe %in% 1:2){
    intervs = c(0, 0.16, 0.32, 0.52, 0.68, 0.88, 1.000001)
  }
  else if (var.to.stripe == 3){
    intervs = (0:6/5.999)^2  
  }
  else if (var.to.stripe == 4){
    intervs = (0:6/5.999)^2
  }
  colormap = mycols[ findInterval(to.plot[, 6], intervs ) ]
  to.plot = data.frame(to.plot, colormap=colormap)

  # Axis labels
  axlab = names(ms.data)[vars.to.plot]

  # Create plot object
  sp3 = scatterplot3d(to.plot[,3],to.plot[,4], to.plot[,5], type="n", 
        grid=FALSE, angle=rot.ang, xlab=axlab[1], ylab=axlab[2], zlab=axlab[3], box=FALSE,
        xlim = c(min(-0.1, min(to.plot[,3])), max(0.1, max(to.plot[,3]))))
  
  # Plot boundary lines
  if(length(setdiff(c(1,2),vars.to.fix)) == 0){
      b1 = to.plot[to.plot$theta < eps,]
      b1 = b1[order(b1$cosphi),]
      b2 = to.plot[to.plot$cosphi == 1,]
      b2 = b2[order(b2$theta),]
      b3 = to.plot[to.plot$theta >= pi/2-eps,]
      b3 = b3[order(b3$cosphi, decreasing = TRUE),]
      b4 = to.plot[to.plot$cosphi == 0,]
      b4 = b4[order(b4$theta, decreasing = TRUE),]
      boundary = rbind(b1,b2,b3,b4)      
      sp3$points3d(boundary[,3], boundary[,4], boundary[,5], type="l", col="black", lwd=2)
    }
    else if(length(setdiff(c(3,4),vars.to.fix))==0){
      b1 = to.plot[to.plot$a3 == 0,]
      b1 = b1[order(b1$a2),]
      b2 = to.plot[to.plot$a2 == 1,]
      b2 = b2[order(b2$a3),]
      b3 = to.plot[ findInterval(to.plot$a2-to.plot$a3, c(-eps,eps)) == 1,]
      b3 = b3[order(b3$a3, decreasing = TRUE),]
      boundary = rbind(b1,b2,b3)      
      sp3$points3d(boundary[,3], boundary[,4], boundary[,5], type="l", col="black", lwd=2)
      if (a.reflect){
        sp3$points3d(boundary[,4], boundary[,3], boundary[,5], type="l", col="black", lwd=2)
      }
    }
  
  # Plot data
  sp3$points3d(to.plot[,3], to.plot[,4], to.plot[,5], type="p", cex=.5, pch=16,
               col = as.character(to.plot$colormap))
  if (a.reflect & length(setdiff(c(3,4),vars.to.fix)) == 0){
    sp3$points3d(to.plot[,4], to.plot[,3], to.plot[,5], type="p", cex=.5, pch=16,
                 col = as.character(to.plot$colormap))
  }
}

plot.fix.3var <- function(ms.data, vars.to.fix, fixed.values, vars.to.plot){
  
  # Select out the data for fixed values
  eps = 1e-8
  to.plot = ms.data[ findInterval( ms.data[,vars.to.fix[1]], 
        fixed.values[1] + c(-eps,eps)) == 1  & findInterval( ms.data[,vars.to.fix[2]], 
        fixed.values[2]+c(-eps,eps))==1 & findInterval( ms.data[,vars.to.fix[3]], 
        fixed.values[3]+c(-eps,eps))==1, c(setdiff(1:4,vars.to.fix), 
                                           vars.to.plot) ]

  # Axis labels
  axlab = names(ms.data)[vars.to.plot]

  # Create plot object
  plot(to.plot[,2],to.plot[,3], type="p", xlab=axlab[1], ylab=axlab[2],
       pch = 16, cex = .6)
  par(mar=c(5,4,4,2))   # I needed to adjust left margin to see the y-axis label
}
  