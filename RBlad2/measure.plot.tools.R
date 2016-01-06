# List of functions in this file:
#
#   plot.mag.ratios()
#   plot.nm.hist()
# 	plot.nm.smooth()
# 	plot.ang.err()
# 	plot.fit.data.naive()
#	plot.fit.data.log()

# Load necessary libraries, scripts and data
library(ggplot2)
library(R.matlab)
source("measure.process.R")

if( !exists("training")){
	training = data.frame(readMat("training.MAT")$training[,,1])
}

# Plots the magnitude ratios of the random and spiral-grid measured training data, as well as a smoothing line.
plot.mag.ratios <- function(){

	# Shape the training data
	training = training[,c("meas","rand.mr","grid.mr")]
	tr.rs = reshape(training,varying=list(c("rand.mr","grid.mr")),direction="long",v.names="mr",times=c("Random","Grid"))

	# Plotting
	g <- ggplot(data=tr.rs, aes(meas,mr))
	g <- g + geom_point() + geom_smooth() + facet_grid(.~time) + xlab("Number of measurements") +
		 ylab("Magnitude ratio") + ggtitle("Ratios of unnormalized measured vectors to actual Bloch vectors")

	return(g)
}

# Histograms for the norms of the measured training data. Shows the overshot measurements using a vertical boundary
plot.nm.hist <- function(){

	# Shape the training data
	if( !exists("train.n")){
		train.n = mm.normalize(training)
	}
	tr.rs = reshape(train.n, varying=list(c("xyz.nm", "rand.nm", "grid.nm")),
		direction="long",v.names="mm.nm",times=c("xyz","Random","Grid"))

	# Plotting
	g <- ggplot(data=tr.rs, aes(mm.nm))
	g <- g + geom_histogram(bins=100) + facet_grid(.~time) + geom_vline(xintercept = 1) + xlab("Magnitude of measured vectors") +
			ggtitle("Magnitudes of measured vectors")

	return(g)
}

# Plots the actual vs measured magnitudes for all three methods, in six different measurement classes. Smoothing lines show how there
# is a hook or levelling off at large measured magnitude. 
plot.nm.smooth <- function(){

	# Load and shape data
	if( !exists("train.n")){
		train.n = mm.normalize(training)
	}
	tr.p = train.n[,c("meas","nmag","xyz.nm")]

	# Bin the measurements and rename
	tr.p = data.frame(tr.p, range = findInterval(train.n$meas,c(100,300,500,1000,1500)))
	tr.p$range = c("Meas < 100", "Meas in (100,300)", "Meas in (300,500)", "Meas in (500,1000)", 
		"Meas in (1000,1500)", "Meas > 1500")[tr.p$range+1]
	tr.p$range = factor(tr.p$range, levels=c("Meas < 100", "Meas in (100,300)", "Meas in (300,500)", 
		"Meas in (500,1000)", "Meas in (1000,1500)", "Meas > 1500"))

	# Plotting
	g <- ggplot(data=tr.p, aes(xyz.nm,nmag))
	g <- g + geom_point() + geom_smooth() + facet_wrap(~range, nrow=2, ncol=3) + 
			xlab("Measured magnitude") + ylab("Actual magnitude") + ggtitle("Actual vs measured magnitudes for the xyz method")

	return(g)
}

# Log-log plots of angular error for four different magnitude levels, according to measurement.
plot.ang.err <- function(){

	# Load and shape data
	if( !exists("train.n")){
		train.n = mm.normalize(training)
	}

	train.n = train.n[,c("meas","xyz.nm", "rand.nm", "grid.nm","xyz.ang.err", "rand.ang.err", "grid.ang.err")]
	ae.rs1 = reshape(train.n, varying=list(c("xyz.ang.err", "rand.ang.err", "grid.ang.err")),
		direction="long",v.names="ang.err",times=c("xyz","Random","Grid"))
	ae.rs2 = reshape(train.n, varying=list(c("xyz.nm", "rand.nm", "grid.nm")),
		direction="long",v.names="mm.nm",times=c("xyz","Random","Grid"))
	ae.rs = merge(ae.rs1,ae.rs2) 
	ae.rs = ae.rs[,c("meas","time","ang.err","mm.nm")]

	# Remove measurements that return the zero vector (angular error is undefined)
	ae.rs = ae.rs[!is.nan(ae.rs$ang.err),]

	# Bin the measured magnitudes
	ae.rs = data.frame(ae.rs, mag.cat = findInterval(ae.rs$mm.nm, (1:3)/4))
	ae.rs.red = aggregate( cbind(ang.err,mm.nm) ~ factor(mag.cat) + factor(meas) + factor(time), data = ae.rs, FUN = mean)
	names(ae.rs.red)[1:3] = c("mag.cat","meas","Method")
	levels(ae.rs.red$mag.cat) <- c("Meas. mag. < .25","Meas. mag. in (.25,.5)","Meas. mag. in (.5,.75)","Meas. mag. > .75")

	# Plotting
	g <- ggplot( data=ae.rs.red, aes(log(as.numeric(meas)), log(ang.err) , color = Method))
	g <- g + geom_point(size=.3, alpha=1) + facet_grid(~mag.cat) + ggtitle("Log-log plot of measurement number vs angular error") +
			ylab("Angular error (log-scale)") + xlab("Measurement number (log-scale)")

	return(g)
}

# Plots the four fit parameters for the segmented fits, for each method, according to measurement.
plot.fit.data.naive <- function(){

	# Load and shape data
	seg.fits = read.table("seg.fits.csv", stringsAsFactors = FALSE)
	seg.fits$method[seg.fits$method=="rand"] = "Random"
	seg.fits$method[seg.fits$method=="grid"] = "Grid"

	seg.fits.long = reshape(seg.fits, varying = list(c("breakpoint", "m1", "m2", "intercept")), direction = "long", v.names = "Fit.parameter",
		times = c("Breakpoint", "Left slope", "Right slope", "Intercept"))

	# Plotting. Turn scales to free, as different parameters have different ranges.
	g <- ggplot(data = seg.fits.long, aes(meas,Fit.parameter)) + geom_point() + facet_grid(time~method, scales="free") +
		ylab("Parameter value") + xlab("Measurements") + ggtitle("Fit parameters for each measurement level")

	return(g)
}

# Plots three of the four fit parameters after transforming them. A few points must be tossed due to undefined logs.
plot.fit.data.log <- function(){

	# Load data
	seg.fits = read.table("seg.fits.csv", stringsAsFactors = FALSE)
	seg.fits$method[seg.fits$method=="rand"] = "Random"
	seg.fits$method[seg.fits$method=="grid"] = "Grid"
	
	# Transform data. Ignore warning about logs of negative values.
	seg.fits.trans = suppressWarnings(data.frame(seg.fits$method, log(seg.fits$meas), log(1-seg.fits$breakpoint), log(seg.fits$intercept), log(1 - seg.fits$m1)) )
	names(seg.fits.trans) = c("method", "logmeas", "log1-bp", "logint","log1-m1")

	# Reshape data frame
	seg.fits.long = reshape(seg.fits.trans, varying = list(c("log1-bp", "logint","log1-m1")), direction = "long", v.names = "Fit.parameter",
		times = c("log(1-breakpoint)", "log(intercept)", "log(1 - left slope)"))

	# Throw out invalid data
	seg.fits.long = seg.fits.long[complete.cases(seg.fits.long) & is.finite(seg.fits.long$Fit.parameter),]

	# Plotting. Add linear fits.
	g <- ggplot(data = seg.fits.long, aes(logmeas,Fit.parameter)) + geom_point() + facet_grid(time~method, scales="free") + 
		geom_smooth(method="lm") + xlab("Measurements (log-scale)") + ylab("Transformed parameter value") +
		ggtitle("Transformed fit parameters for each measurement level")

	return(g)
}