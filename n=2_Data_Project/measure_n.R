# This function simulates a measurement of n given its actual value
# The number of trials is given as an input variable

measure_n <- function(no_tri, actual_n){
 prob = (actual_n + 1)/2        # find probabilities
 return( 2*rbinom(3, no_tri, prob)/no_tri - 1) # generate measurement
}