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

measure.test <- function(no.n, no.meas, mag.range = 0:1 ){

  # Compute size of data frame and intialize  
  nr = no.n*length(no.meas)
  results = data.frame( meas = double(nr),
            n1 = double(nr), n2 = double(nr), n3 = double(nr), nmag = double(nr),
            xyz.n1 = double(nr), xyz.n2 = double(nr), xyz.n3 = double(nr), xyz.nm = double(nr),
            xyz.mag.err = double(nr), xyz.ang.err = double(nr), xyz.tot.err = double(nr),
            rand.n1 = double(nr), rand.n2 = double(nr), rand.n3 = double(nr), 
            rand.nm = double(nr), rand.mr = double(nr), rand.ang.err = double(nr), 
            grid.n1 = double(nr), grid.n2 = double(nr), grid.n3 = double(nr),
            grid.nm = double(nr), grid.mr = double(nr), grid.ang.err = double(nr))
  
  # Loop over measurement numbers
  for (l in 1:length(no.meas)){
    mm = no.meas[l]

    # Loop over vectors
    for (j in 1:no.n){
      print(paste(mm, " measurements, vector ", j, sep=""))

      # Generate actual Bloch vector
      n = rand.n()
      n.mag = mag(n)*(mag.range[2]-mag.range[1])+mag.range[1]
      n = n*n.mag/mag(n)
      
      # Measure using the x-y-z method
      n.xyz = c( measure(round(mm/3), n, c(1,0,0)),
                 measure(round(mm/3), n, c(0,1,0)),
                 measure(mm - 2*round(mm/3), n, c(0,0,1)))
      n.xyz.mag = mag(n.xyz)

      # Calculate errors in the x-y-z method
      n.xyz.me = n.xyz.mag - n.mag
      n.xyz.ae = 1-(n %*% n.xyz)/ n.mag / n.xyz.mag
      n.xyz.te = mag(n-n.xyz)
      
      # Measure using the random-axis method
      n.rand = c(0,0,0)
      for (k in 1:mm){
        m = rand.n()
        m = m/mag(m)
        n.rand = n.rand + measure(1, n, m)*m
      }
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
        n.grid = n.grid + measure(1,n,m)*m
      }
      n.grid.mag = mag(n.grid)
      
      # Compute angular error
      n.grid.ae = 1 - (n %*% n.grid)/ n.mag / n.grid.mag

      # Record results      
      results[(l-1)*no.n + j, ] = c(mm, n, n.mag, n.xyz, n.xyz.mag, n.xyz.me, n.xyz.ae, n.xyz.te, 
                n.rand, n.rand.mag, n.rand.mag/n.mag, n.rand.ae, 
                n.grid, n.grid.mag, n.grid.mag/n.mag, n.grid.ae)
    }
  }    
  return(results)
}

# This function processes the raw data for Bloch ball measurements. The goals are:
#   (1) model and compare the angular error for the three methods as it varies with measurement number
#   (2) model the magnitude ratios for the random-vector and spiral-grid methods
#   (3) using the aforementioned magnitude ratios, calculate estimates for the random-vector and spiral-grid methods
#   (4) model and compare the magnitude error for the three methods as it varies with measurement number

mm.proc <- function(raw.df){

  # Initialize summary data frame
  raw.df = raw.df[order(raw.df$meas),]
  nsr = length(unique(raw.df$meas))
  summ.df = data.frame(meas = sort(unique(raw.df$meas)), 
          av.xyz.ang.err = double(nsr), av.rand.ang.err = double(nsr), av.grid.ang.err = double(nsr), 
          sd.xyz.ang.err = double(nsr), sd.rand.ang.err = double(nsr), sd.grid.ang.err = double(nsr), 
          av.rand.mag.rat = double(nsr), av.grid.mag.rat = double(nsr), 
          sd.rand.mag.rat = double(nsr), sd.grid.mag.rat = double(nsr), 
          av.xyz.mag.err = double(nsr), av.rand.mag.err = double(nsr), av.grid.mag.err = double(nsr), 
          sd.xyz.mag.err = double(nsr), sd.rand.mag.err = double(nsr), sd.grid.mag.err = double(nsr), 
          av.xyz.abmag.err = double(nsr), av.rand.abmag.err = double(nsr), av.grid.abmag.err = double(nsr), 
          sd.xyz.abmag.err = double(nsr), sd.rand.abmag.err = double(nsr), sd.grid.abmag.err = double(nsr), 
          av.xyz.tot.err = double(nsr), av.rand.tot.err = double(nsr), av.grid.tot.err = double(nsr), 
          sd.xyz.tot.err = double(nsr), sd.rand.tot.err = double(nsr), sd.grid.tot.err = double(nsr))
          
  # Summarize angular errors and magnitude ratios
  summ.df[,c(2:4,8,9)] = aggregate(raw.df[,c(11,18,24,17,23)], by = list(raw.df[,1]), mean )[,2:6] 
  summ.df[,c(5:7,10,11)] = aggregate(raw.df[,c(11,18,24,17,23)], by = list(raw.df[,1]), sd )[,2:6]
  
  # Adjust random-vector and spiral-grid estimates
  
    rsl = rep(summ.df[,8],each=10000)
    gsl = rep(summ.df[,8],each=10000)
    raw.df[,13:16] = raw.df[,13:16]/rsl
    raw.df[,19:22] = raw.df[,19:22]/gsl
    raw.df = raw.df[,c(1:16,18:22,24)]  # Discard magnitude ratios
    
    # Calculate magnitude and total errors for the random-vector and spiral-grid methods
    rme = raw.df[,16] - raw.df[,5]
    gme = raw.df[,21] - raw.df[,5]
    rte = mag(raw.df[,13:15] - raw.df[2:4])
    gte = mag(raw.df[,18:20] - raw.df[2:4])

    raw.df = data.frame(raw.df[,1:16], rand.mag.err = rme, rand.tot.err = rte, 
                        raw.df[,17:21], grid.mag.err = gme, grid.tot.err = gte)

    summ.df[,c(12:14,24:26)] = aggregate(raw.df[,c(10,17,24,12,18,25)], by = list(raw.df[,1]), mean )[,2:7] 
    summ.df[,c(15:17,27:29)] = aggregate(raw.df[,c(10,17,24,12,18,25)], by = list(raw.df[,1]), sd )[,2:7]
    summ.df[,18:20] = aggregate( abs(raw.df[,c(10,17,24)]), by = list(raw.df[,1]), mean )[,2:4] 
    summ.df[,21:23] = aggregate( abs(raw.df[,c(10,17,24)]), by = list(raw.df[,1]), sd )[,2:4] 
    
    return( list( Summarized = summ.df, Processed = raw.df) ) 
}


mag.err.dep <- function(raw.df, incr=.1){
  
  x = unique(raw.df$meas)
  methods = c("xyz", "rand", "grid")

  Ly = ceiling(max(cbind(raw.df$xyz.nm, raw.df$rand.nm, raw.df$grid.nm)) / incr)  
  dep = data.frame( meas = rep(x, each = Ly), incr = rep(incr*(1:Ly), length(x)), 
        av.xyz.mag.err = double(length(x)*Ly),
        av.rand.mag.err = double(length(x)*Ly),
        av.grid.mag.err = double(length(x)*Ly))
  
  for (k in 1:3){
    zp = raw.df[,paste(methods[k],"mag.err",sep=".")]

    for (j in 1: (length(x)*Ly)){
        z = zp[ raw.df$meas == x[ ceiling(j/Ly)] & 
            findInterval(raw.df[,paste(methods[k],"nm",sep=".")], 
                    incr*( j%%Ly + -1:0)) == 1 ] 
        dep[j,2+k] = mean(z)
    }
  }
  return(dep)
}

two.line.fit <- function(raw.df){
  
  methods = c("xyz","rand","grid")
  meas = sort(unique(raw.df$meas))
  init.guess = c(1,0,0,1)

  nr = 3*length(meas)  
  fit.df = data.frame(meth = rep(methods, each = length(meas)), meas = rep(meas,3), 
            xc = double(nr), yc = double(nr), 
            m1 = double(nr), m2 = double(nr))

  for (j in 1:nr){
    
    cm = meas[ ((j-1) %% length(meas)) +1]
    cmth = methods[ ceiling(j/length(meas)) ]
    print( paste("Fitting data for", cm, "measurements and method", cmth, sep = " "))

    pts = raw.df[ raw.df$meas == cm, ]
    me.sel = names(pts) == paste("av",cmth,"mag","err",sep=".")
  
    pts = data.frame(incr = pts$incr, mag.err=pts[,me.sel])
    pts = pts[!is.nan(pts$mag.err),]
    pts = pts[pts$incr > .2, ]
    pp = optim(init.guess, function(params){return(tlf.cost.fn(params, pts))})$par  
    
    fit.df[j, 3:6] = pp
  }
  return(fit.df)
}

mms.adjust <- function(raw.df, gtlf, method="xyz"){

  meas = unique(raw.df$meas)
  nr = dim(raw.df)[1]
  ret.df = data.frame(m.n1=double(nr),m.n2=double(nr),m.n3=double(nr),m.nm=double(nr),
                      m.me=double(nr),m.ae=double(nr),m.te=double(nr))
  count = 1
  for (j in 1:length(meas)){
    w = meas[j]
    m1 = gtlf$m1.coeff*w^gtlf$m1.pow
    m2 = gtlf$m2
    if (method =="xyz"){
      xc = 1
      yc = gtlf$yc.coeff*w^gtlf$yc.pow
    }
    else{
      xc = 1 - gtlf$xc.coeff*(w-gtlf$xc.shift)^gtlf$xc.pow      
      yc = 0
    }
    spec.df = raw.df[raw.df$meas == w,]

    for (k in 1:dim(spec.df)[1]){
      
      x = spec.df[k,paste(method,"nm",sep=".")]
      if (x < xc){
        y =  m1*(x-xc) + yc
      }
      else{
        y = m2*(x-xc) + yc
      }
      mult = 1/(1 + y)

      rn = row.names(spec.df[k,])
      for (l in 1:3){
        ret.df[count,l] = raw.df[rn,paste(method,".n",l,sep="")]*mult
      }
      ret.df[count,4] = mag(ret.df[count,1:3])
      ret.df[count,5] = ret.df[count,4] - raw.df[rn,"nmag"] 
      ret.df[count,6] = raw.df[rn,paste(method,"ang.err",sep=".")]
      ret.df[count,7] = mag(ret.df[count,1:3]-raw.df[rn,c("n1","n2","n3")])
      count = count + 1
    }
  }
  return(ret.df)
}

tlf.cost.fn <- function(params, pts){
  
  pts = pts[!is.nan(pts[,2]),]
  
  cost = 0
  xc = params[1]
  yc = params[2]
  m1 = params[3]
  m2 = params[4]

  for (j in 1:dim(pts)[1]){
    x = pts[j,1]
    y = pts[j,2]
    if (x < xc){
      cost = cost + (y - m1*(x-xc) - yc )^2
    }
    else {
      cost = cost + (y - m2*(x-xc) - yc)^2
    }
  }
  
  return(cost)
}

global.two.line.fit <- function(raw.df, method="xyz"){
  
  if (method=="xyz"){
    init.guess = c(2, -1, 2, -.75, .85)
  }
  else if (method=="rand"){
    init.guess = c(3,-1, 0, 2,-.75,.75)
  }
  else if (method=="grid"){
    init.guess = c(3,-1, 0, 2,-.75,.75)
  }
  
  fit.df = NULL
  
  me.sel = names(raw.df) == paste("av",method,"mag","err",sep=".")

  pts = data.frame(meas=raw.df$meas, incr = raw.df$incr, mag.err=raw.df[,me.sel])
  pts = pts[pts$incr > 0.2,]
  pp = optim(init.guess, function(params){return(global.tlf.cost.fn(params, pts, 
                                                      method))})$par  

  if (method == "xyz"){
    fit.df$yc.coeff = pp[1]
    fit.df$yc.pow = pp[2]
    fit.df$m1.coeff = pp[3]
    fit.df$m1.pow = pp[4]
    fit.df$m2 = pp[5]
  }    
  else {
    fit.df$xc.coeff = pp[1]
    fit.df$xc.pow = pp[2]
    fit.df$xc.shift = pp[3]
    fit.df$m1.coeff = pp[4]
    fit.df$m1.pow = pp[5]
    fit.df$m2 = pp[6]
  }    

  
  return(fit.df)
}

global.tlf.cost.fn <- function(params, pts, method="xyz"){
  pts = pts[!is.nan(pts[,3]),]

  cost = 0
  coord.coeff = params[1]
  coord.pow = params[2]
  
  if (method=="xyz"){
    m1.coeff = params[3]
    m1.pow = params[4]
    m2 = params[5]
  }
  else{
    coord.shift = params[3]
    m1.coeff = params[4]
    m1.pow = params[5]
    m2 = params[6]
  }

  for (j in 1:dim(pts)[1]){
    x = pts[j,2]
    y = pts[j,3]
    if(method=="xyz"){
      xc = 1  
      yc = coord.coeff*(pts[j,1])^coord.pow
    }
    else{
      xc = 1 - coord.coeff*(pts[j,1]-coord.shift)^coord.pow
      yc = 0
    }
    m1 = m1.coeff*pts[j,1]^m1.pow
    if (x < 1){
      cost = cost + pts[j,1]*( y - m1*(x-xc) - yc )^2
    }
    else {
      cost = cost + pts[j,1]*( y - m2*(x-xc) - yc )^2
    }
  }
  return(cost)
}


# This function simulates a measurement of n given its actual value
# The number of trials is given as an input variable. The measurements are made along
# the x,y,z axes.

measure <- function(no.meas, actual.n, direction){
  direction = direction / sqrt(sum(direction^2))
  prob = (actual.n %*% direction + 1)/2        # find probabilities
  return( 2*rbinom(1, no.meas, prob)/no.meas - 1) # generate measurement
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

# Apparently, R has no function for computing the magnitude of a vector (!?), 
# so here it is. Also works

mag <- function(n){
  if (class(n) == "data.frame"){
    return(sqrt(rowSums(n^2)))
  }
  else{
    return(sqrt(sum(n^2)))
  }
}
