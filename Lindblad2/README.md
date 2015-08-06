# Lindblad2

This is a program for analyzing the fast control of a 2d Lindblad quantum control system. The goals are: (1) Find two "threads" in the Bloch ball: the points on concentric surfaces where the (outward) radial velocity is maximized or minimized. The minimum radial velocity is always non-positive; the maximum velocity is positive at the center, and decreases monotonically as one moves outward. At some point before the outermost sphere, the maximum velocity crosses zero - this point is called the "horizon". (2) Find the "chimney", the region where outward velocity is positive.  

In some special cases, there are bifurcations (either thread may split into three, or a third thread may pop up independently, or a plane of CPs emerges from a thread). These threads are detected by the program and dealt with accordingly. If the maximum-velocity thread bifurcates, there may be an alternate horizon (actually, a pair of equal-radius horizons). 

### Installation
 
This program depends on Armadillo, a C++ linear algebra package:

http://arma.sourceforge.net/

Conrad Sanderson. 
*Armadillo: An Open Source C++ Linear Algebra Library for Fast Prototyping and Computationally Intensive Experiments.* 
Technical Report, NICTA, 2010.

Armadillo depends in turn on LAPACK and BLAS.

If all three are properly installed, simply run `make` from the *src* directory. Additionally the user should have gnuplot installed in order to process the program results. The .gnu file that is produced assumes you have the "windows" terminal. If you want the "wxt" or different, the files must be altered. Additionally, the minimum threads are dashed but may not appear so if the sampling frequency is too high. To get dashed threads, the .gnu files must be altered.

### Usage

Lindblad2 is an interactive program from the terminal. There are five options for producing a system: 

1. The program can choose a random system (this is the appropriate choice for users unfamiliar with the math).
2. The six system parameters can be stored as a file.
3. The six system parameters can be input at the terminal.
4. The nine parameters of the GKS matrix (using the Pauli matrix basis) can be input at the terminal.
5. Similar to 4, but the GKS matrix uses the +/-/z matrix basis.

The user will be prompted for a discretization number. The recommended number is 1000. If the thread curvature is high (which sometimes happens near the center), the user can increase to 10,000 or 100,000, but in the latter case the run-time can be a few minutes.

Additionally, the user must give a file name base which the program uses to name the output files.

### Output

For each system, the program produces ten files: 

* `***_n.dat` contains the thread co-ordinates. Bifurcated threads are also stored there, but the formatting is not important as the user is not meant to deal directly with this file. 
* `***_w.dat` contains the chimney co-ordinates. The algorithm does not deal well with the top of the chimney, and so the terminal point is inserted manually as the horizon on the maximizing thread, only if the lines are sufficiently close. Occasionally however, this is not the right thing to do! Some chimney lines terminate at a different location, usually if the thread has high curvature. In this case, the lines don't terminate and there is a hole in the chimney.
* `***_n_plot.gnu` processes the thread .dat file for gnuplot. It uses the wxt terminal to plot. The line color and types can vary between machines, so the user will have to alter this file if their machine does not produce satisfactory results. All threads are meant to be in black; the maximal thread should be solid, the minimal thread dotted, and any bifurcated threads dashed. If the discretization number is much higher than 1000, then the dotted/dashed lines may be too "dense" to look dotted or dashed, in which case the user must tweak the .gnu file.
* `***_f.dat` contains the radial velocities as functions of r. 
* `***_f_plot.gnu` processes the above file. Main threads are solid curves, and bifurcated threads are dotted.
* `***_e.dat` contains the estimated errors for the threads (not for bifurcated threads as these cases can be solved exactly).
* `***_e_plot.gnu` processes the above file with gnuplot.
* `***_we.dat` contains the estimated errors for the chimney.
* `***_we_plot.gnu` processes the above file with gnuplot.
* The six system parameters are output in `***_param.txt`, unless the system parameters were input as a file.
