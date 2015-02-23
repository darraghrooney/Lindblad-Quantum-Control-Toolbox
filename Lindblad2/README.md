# Lindblad2

This is a program for analyzing the fast control of a 2d Lindblad quantum control system. The goal is to essentially to find two "threads" in the Bloch ball: the points on concentric surfaces where the (outward) radial velocity is maximized or minimized. The minimum radial velocity is always non-positive; the maximum velocity is positive at the center, and decreases monotonically as one moves outward. At some point before the outermost sphere, the maximum velocity crosses zero - this point (the "horizon") is marked with a large black dot. 

In some special cases, there are bifurcations (either thread may split into three, or a third thread may pop up independently). These threads are detected by the program and dealt with. If the maximum-velocity thread bifurcates, there will be an alternate horizon (actually, a pair of equal-radius horizons). 

## Installation:
 
This program depends on Armadillo, the C++ linear algebra package:

http://arma.sourceforge.net/

Conrad Sanderson. 
*Armadillo: An Open Source C++ Linear Algebra Library for Fast Prototyping and Computationally Intensive Experiments.* 
Technical Report, NICTA, 2010.

Armadillo depends in turn on LAPACK and BLAS.

If all three are properly installed, simply run `make` from the *src* directory. Additionally the user should have gnuplot installed in order to process the program results.

## Usage

Lindlad2 is an interactive program from the terminal. There are five options for producing a system: 
1. The program can choose a random system (this is the appropriate choice for users unfamiliar with the math).
2. The six system parameters can be stored as a file.
3. The six system parameters can be input at the terminal.
4. The nine parameters of the GKS matrix (using the Pauli matrix basis) can be input at the terminal.
5. Similar to 4, but the GKS matrix uses the +/-/z matrix basis.

The user will be prompted for a discretization number. The recommended number is 1000. If the thread curvature is high (which sometimes happens near the center), the user can increase to 10,000 or 100,000, but in the latter case the run-time can be a few minutes.

Additionally, the user must give a file name base which the program uses to name the output files.

## Output

For each system, the program produces seven files: 
* `***_n.dat` contains the thread co-ordinates. Horizon co-ordinates as well as bifurcated threads are also stored there, but the formatting is not important as the user is not meant to deal. 
* `***_n_plot.gnu` processes the thread .dat file for gnuplot. It uses the wxt terminal to plot. The line color and types can vary between machines, so the user will have to alter this file if their machine does not produce satisfactory results. All threads are meant to be in black; the maximal thread should be solid, the minimal thread dotted, and any bifurcated threads dashed.
* `***_f.dat` contains the radial velocities as functions of r. 
* `***_f_plot.gnu` processes the above file. Main threads are solid curves, and bifurcated threads are dotted.
* `***_e.dat` contains the estimated errors for the threads (not for bifurcated threads as these cases can be solved exactly).
* `***_e_plot.gnu` processes the above file with gnuplot.
* The six system parameters are output in `***_param.txt`, but not if the system parameters were input as a file.