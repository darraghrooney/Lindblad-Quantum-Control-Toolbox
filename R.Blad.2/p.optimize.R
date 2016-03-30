library(R.matlab)
train.set = data.frame(readMat("assess.train.MAT")$assess[,,1])

library(deldir)
source("assessment.R")

p.optimize <- function(){

  train.size = dim(train.set)[1]
  tol = -2

  nan.pts = is.nan(train.set$trace.metric) & is.nan(train.set$diag.metric)
  bad.pts1 = log(1-train.set$trace.metric) > tol
  bad.pts1[is.na(bad.pts1)] = TRUE
  bad.pts2 = log(1-train.set$diag.metric) > tol
  bad.pts3 = log(1-train.set$off.diag.metric) > tol
  bad.pts4 = log(1-train.set$b.metric) > tol
  bad.set = train.set[nan.pts | bad.pts1 | bad.pts2 | bad.pts3 | bad.pts4,]
  bad.size = dim(bad.set)[1]
  
  dd = deldir(bad.set[,c("p1", "p2")])
  plot.tile.list(tile.list(dd),border = 2,xlim=c(1,3),ylim=c(0,.5))
  p.opt = best.pt(dd)
  
  return(list(bad.set = bad.set, p.opt = p.opt))
}

best.pt <- function(dd){
  deltri = triang.list(dd)
  best.pt = c(0,0,0)
  for (j in 1:length(deltri)){
    new.pt = circumradius(deltri[[j]])
    if (new.pt$cr > best.pt[3] & new.pt$cc[1] > 1  & new.pt$cc[1] < 3
        & new.pt$cc[2] > 0 & new.pt$cc[2] < 0.5){
      best.pt = c(new.pt$cc, new.pt$cr)
    }
  }
  return(best.pt)
}

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
