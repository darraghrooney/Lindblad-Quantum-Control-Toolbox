# This function generates (pieces of) the "muscape" ... the Lagrangian multipliers
# for optimizing the horizon for n=2 dissipative quantum systems. The original parameter
# space is 9D, but after accounting for various symmetries, it can be reduced to
# 4D. Two of these are in the 3-vector a. The largest component of a is normalized
# to 100, and the other two are the parameters 100 >= a2 >= a3 >= 0. The other two
# come from the 3-vector b. mu is invariant under b->c*b, which leaves two independent
# variables: cos(phi) and theta, where b1=sqrt(a2a3)*sin(theta)*cos(phi), 
# b2 = sqrt(a1a3)*sin(phi)*sin(theta) and b3 = sqrt(a1a2)*cos(phi).

# Input values are a_disc and ang_disc which discretize the ranges a2,a3 in [0,100] for
# a_disc, and cp in [0,1] and theta in [0,pi/2] for ang_disc. The entire muscape will
# be generated in pieces, so we allow a2.inds and a3.inds to specify which values are
# done. Note that a2.inds must be in (0,100]/a_disc and a3.inds in [0,a2.inds].

# The resulting data frame contains the parameters, the mu, the corresponding 
# horizon co-ordinates, and the horizon magnitude

muscape_generate <- function(a_disc,ang_disc,a2.inds,a3.inds){
  
  # Figure out how many rows will be in the data frame
  ac = 0
  for (j in a2.inds){
    for (j in a3.inds){
      ac = ac+1
    }
  }
  no_pts = ang_disc^2*ac
  
  # Initialize data frame
  muscape = data.frame(a2=double(no_pts),a3=double(no_pts),cp=double(no_pts),
         th=double(no_pts),mu=double(no_pts),nM1=double(no_pts),nM2=double(no_pts),
         nM3=double(no_pts),nMmag=double(no_pts))
  count = 1

  # We need the functions horizon_find() and b_from_angles()
  source("horizon_find.R")
  
  # Sweep over a2, a3 indices
  for (j1 in a2.inds){
    for (j2 in a3.inds){
      a=c(100,100/a_disc*j1,100/a_disc*j2)
      
      # Sweep over angles
      for (j3 in 0:ang_disc){
        for (j4 in 0:ang_disc){

          # Calculate data
          cp = j3/ang_disc
          th = j4/ang_disc/2*pi
          b= b_from_angles(a,cp,th,1)
          mus = horizon_find(a,b)
          muscape[count,] = c(a[2],a[3],cp,th,mus)
          count = count+1
        }
      }
    }
  }
  return(muscape)
}


# This function takes (pieces of) the original muscape, and calculates the necessary
# Hamiltonian for each horizon. These Hamiltonians are added to the data frame.

muscape_add_hams <- function(msc){

  # Determine size of muscape piece
  nrows=dim(msc)[1]
  if (dim(msc)[2] == 9){
    
    # Add necessary columns if not present
    msc = data.frame(msc, hM1=double(nrows), hM2=double(nrows), hM3=double(nrows))
  }
  
  # Sweep over rows
  for (j in 1:nrows){
    
    # Calculate necessary data and insert Hamiltonian into data frame
    msrow = as.matrix(msc[j,])
    a = as.vector(c(100, msrow[1], msrow[2]))
    b = b_from_angles(a, msrow[3],msrow[4],1)
    msc[j, 10:12] = h_from_n(msrow[6:8], a, b)
  }
  return(msc)
}