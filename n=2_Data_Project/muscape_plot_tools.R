# These are functions to produce various plots for the project

# Plots the discretization of a-space.
plot_as <- function(ms_data){

  # Create data frame
  a_data = data.frame(unique(ms_data[,1:2]))
  nr = dim(a_data)[1]
  a_data = data.frame(a_data, color = character(nr))

  # Assign colors by striping according to the a3 variable
  mycols = c("magenta", "blue","green", "red", "darkorange", "violet")
  a_data$color = mycols[findInterval(a_data$a3, c(0,20,40,60,80,100))]

  # Plot
  plot(c(0,100,100,0), c(0,0,100,0), xlab="a2", ylab="a3", type="l")  # Boundary
  points(a_data$a2,a_data$a3, col = a_data$color, pch = 19, cex = .5) # Data
  par(mar=c(5,4,4,2))   # I needed to adjust left margin to see the y-axis label
}

# Plots the discretization of angle-space
plot_angles <- function(ms_data, rot_ang = 70){
  
  # Create data frame
  angle_data = data.frame(unique(ms_data[,3:4]))
  nr = dim(angle_data)[1]
  sp = sqrt(1-angle_data$cp^2)
  angle_data = data.frame(angle_data, bx = sp*cos(angle_data$theta), 
                  by = sp*sin(angle_data$theta), color = character(nr))
  
  # Assign colors according to cosphi variable
  mycols = c("magenta", "blue","green", "red", "darkorange", "violet")
  angle_data$color = mycols[findInterval(angle_data$cp, 0:5/4.999)]
  
  # Plot
  library(scatterplot3d)
  ad = 0:99/200*pi
  sp3 = scatterplot3d(angle_data$bx,angle_data$by, angle_data$cp,  
            grid=FALSE, box = FALSE, type="n", angle = rot_ang, xlab = "b1",
            ylab="b2", zlab="b3")
  # Boundary
  sp3$points(c(cos(ad),rep(0,99),sin(ad)), c(sin(ad),cos(ad),rep(0,99)), 
                      c(rep(0,99),sin(ad),cos(ad)), type="l")
  # Data
  sp3$points(angle_data$bx,angle_data$by, angle_data$cp, col = angle_data$color, 
         pch = 19, cex = .5)
}

# This functions plots three variables (out of a possible twelve) while
# fixing either a2 and a3, or cosphi and theta, and
# sweeping over the complementary pair. vars_to_fix should be either 
# c(1,2) or c(3,4), and vars_to_plot should be a triple from 1:12
# (that is disjoint from vars_to_fix). The plots are striped, so 
# var_to_stripe should be in 1:4 but not in vars_to_fix

plot_fix_2var <- function(ms_data, vars_to_fix, fixed_values, vars_to_plot,
    var_to_stripe, rot_ang=150){
  
  library(scatterplot3d)
  
  # Select out the data for fixed values
  eps = 1e-8
  to_plot = ms_data[ findInterval( ms_data[,vars_to_fix[1]], 
        fixed_values[1] + c(-eps,eps)) == 1  &
        findInterval( ms_data[,vars_to_fix[2]], fixed_values[2]+c(-eps,eps))==1, 
        c(setdiff(1:4,vars_to_fix), vars_to_plot,var_to_stripe) ]

  # Assign colors
  colormap = character(dim(to_plot)[1])
  mycols = c("magenta", "blue","green", "red", "darkorange", "violet")
  if (var_to_stripe %in% 1:2){
    intervs = c(0, 16, 32, 52, 68, 88,100.0001)
  }
  else if (var_to_stripe == 3){
    intervs = 0:6/5.999  
  }
  else if (var_to_stripe == 4){
    intervs = (0:6/5.999)*pi/2
  }
  colormap = mycols[ findInterval(to_plot[, 6], intervs ) ]
  to_plot = data.frame(to_plot, colormap=colormap)

  # Axis labels
  axlab = names(ms_data)[vars_to_plot]
  axlab[axlab == "cp"] = "cosphi"
   
  # Create plot object
  sp3 = scatterplot3d(to_plot[,3],to_plot[,4], to_plot[,5], type="n", 
        grid=FALSE, angle=rot_ang, xlab=axlab[1], ylab=axlab[2], zlab=axlab[3], box=FALSE,
        xlim = c(min(-0.1, min(to_plot[,3])), max(0.1, max(to_plot[,3]))))
  
  # Plot boundary lines
  if(length(setdiff(c(1,2),vars_to_fix)) == 0){
      b1 = to_plot[to_plot$theta == 0,]
      b1 = b1[order(b1$cp),]
      b2 = to_plot[to_plot$cp == 1,]
      b2 = b2[order(b2$theta),]
      b3 = to_plot[to_plot$theta >= pi/2-eps,]
      b3 = b3[order(b3$cp, decreasing = TRUE),]
      b4 = to_plot[to_plot$cp == 0,]
      b4 = b4[order(b4$theta, decreasing = TRUE),]
      boundary = rbind(b1,b2,b3,b4)      
      sp3$points3d(boundary[,3], boundary[,4], boundary[,5], type="l", col="black", lwd=2)
    }
    else if(length(setdiff(c(3,4),vars_to_fix))==0){
      b1 = to_plot[to_plot$a3 == 0,]
      b1 = b1[order(b1$a2),]
      b2 = to_plot[to_plot$a2 == 100,]
      b2 = b2[order(b2$a3),]
      b3 = to_plot[ findInterval(to_plot$a2-to_plot$a3, c(-eps,eps)) == 1,]
      b3 = b3[order(b3$a3, decreasing = TRUE),]
      boundary = rbind(b1,b2,b3)      
      sp3$points3d(boundary[,3], boundary[,4], boundary[,5], type="l", col="black", lwd=2)
    }
  
  # Plot data
  sp3$points3d(to_plot[,3], to_plot[,4], to_plot[,5], type="p", cex=.5, pch=16,
               col = as.character(to_plot$colormap))
}

plot_fix_3var <- function(ms_data, vars_to_fix, fixed_values, vars_to_plot){
  
  # Select out the data for fixed values
  eps = 1e-8
  to_plot = ms_data[ findInterval( ms_data[,vars_to_fix[1]], 
        fixed_values[1] + c(-eps,eps)) == 1  & findInterval( ms_data[,vars_to_fix[2]], 
        fixed_values[2]+c(-eps,eps))==1 & findInterval( ms_data[,vars_to_fix[3]], 
        fixed_values[3]+c(-eps,eps))==1, c(setdiff(1:4,vars_to_fix), 
                                           vars_to_plot) ]

  # Axis labels
  axlab = names(ms_data)[vars_to_plot]
  axlab[axlab == "cp"] = "cosphi"
  
  # Create plot object
  plot(to_plot[,2],to_plot[,3], type="p", xlab=axlab[1], ylab=axlab[2],
       pch = 16, cex = .6)
  par(mar=c(5,4,4,2))   # I needed to adjust left margin to see the y-axis label
}
  