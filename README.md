# Lindblad Quantum Control Toolbox

This repository contains C++ and Octave/MATLAB code for analyzing fast quantum control systems subject to Lindblad dissipation. Right now, these tools consist of:

* Lindblad2, an interactive C++-based program for finding the maximizing and minimizing threads in the Bloch ball
* MATLAB/Octave code for analyzing systems on two-dimensional Hilbert spaces. Primarily for optimizing control at each radius. (Also there is some code for sweeping over the parameter space.) The optimization code here is different from Lindblad2 in that each radius is optimized independently, whereas Lindblad2 uses the known optimization at r=0, then uses an ODE to track the optimization outwards.
* MATLAB/Octave code for analyzing systems on three-dimensional Hilbert spaces. 
* MATLAB/Octave code for analyzing systems on four-dimensional Hilbert spaces
* MATLAB/Octave code for calculating 

Details about tools are described in the corresponding README's in each directory.
