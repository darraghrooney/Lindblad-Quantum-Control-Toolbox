/* 
 * TwinBlad
 *
 * Copyright (c) 2014 - 2015 Patrick Rooney (darraghrooney@gmail.com)
 *
 * This file is part of TwinBlad.
 *
 * TwinBlad is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TwinBlad is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TwinBlad.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "LSystem2d.h"
#define CHIMNEY_TOL 0.003 
#define CHIMNEY_LINES 36
#define SPOKES 48

using namespace std;
using namespace arma;

// Default constructor assumes there is one L-op and it is equal to sigma_+
LSystem2d::LSystem2d(void)
{	avec = vec("1 1 0");
	bvec = vec("0 0 2");
	rotation = eye(3,3);}

// Constructor for avec and bvec as input
LSystem2d::LSystem2d(vec * inavec, vec * inbvec)
{		
	// Test that the a_j's are ordered and non-negative
	if ((*inavec)(0) < (*inavec)(1)-1E-15 || 
		(*inavec)(1) < (*inavec)(2)-1E-15 || 
		(*inavec)(2) < -1E-15)
		throw 1; 

	// Test that the bvec is small enough	
	else if (4*(*inavec)(0)*(*inavec)(1)*(*inavec)(2) <= 
			(*inavec)(1)*(*inbvec)(1)*(*inbvec)(1)+
			(*inavec)(2)*(*inbvec)(2)*(*inbvec)(2)+
			(*inavec)(0)*(*inbvec)(0)*(*inbvec)(0) -1E-15)
		throw 2;

	// Test that the bvec is non-zero
	else if (norm(*inbvec, 2) == 0)
		throw 3;

	// Assign input to member variables. Default rotation matrix is identity
	else {	
		avec = *inavec; 
		bvec = *inbvec; 
		rotation = eye(3,3);
	}
}

// Constructor for input as avec, bvec and rotation matrix
LSystem2d::LSystem2d(vec * inavec, vec * inbvec, mat * inrot)
{		
	// Test that a_j's are ordered and non-negative
	if (	(*inavec)(0) < (*inavec)(1)-1E-15 || 
		(*inavec)(1) < (*inavec)(2)-1E-15 || 
		(*inavec)(2) < -1E-15)
		throw 1;

	// Test that bvec is small enough	
	else if (4*(*inavec)(0)*(*inavec)(1)*(*inavec)(2) <= 
			(*inavec)(1)*(*inbvec)(1)*(*inbvec)(1)+
			(*inavec)(2)*(*inbvec)(2)*(*inbvec)(2)+
			(*inavec)(0)*(*inbvec)(0)*(*inbvec)(0) -1E-15)
		throw 2;

	// Test that bvec is non-zero
	else if (norm(*inbvec, 2) == 0)
		throw 3;

	// Assign input to member variables	
	else {	
		avec = *inavec; 
		bvec = *inbvec; 
		rotation = *inrot;
	}
}

// Constructor for random Lindblad system
LSystem2d::LSystem2d(char isr)
{	
	srand ((unsigned int)time(NULL));
	if (isr == 'r'){
		
		// a_1 is normalized to 100
		avec = vec("100 0 0");

		// a_2 is random on [0,100]
		avec(1)=100*(double)rand()/(double)RAND_MAX;

		// a_3 is random on [0,a_2]
		avec(2)=avec(1)*(double)rand()/(double)RAND_MAX;

		// Calculate random angles on the sphere
		double theta = 2*PI*(double)rand()/(double)RAND_MAX;
		double phi = PI*(double)rand()/(double)RAND_MAX;

		// Calculate random quasi-norm of bvec
		double beta = (double)rand()/(double)RAND_MAX;

		// Calculate bvec
		bvec = zeros<vec>(3);
		bvec(0) = 2*sqrt(avec(1)*avec(2))*sin(phi)*cos(theta);
		bvec(1) = 2*sqrt(avec(0)*avec(2))*sin(phi)*sin(theta);
		bvec(2) = 2*sqrt(avec(1)*avec(0))*cos(phi); 
		bvec *= beta;

		// Rotation is ignored
		rotation = eye(3,3);
	} 
	else
		// No other character input makes sense
		throw 4;
}

// Accessor functions for the member variables (none for the rotation matrix)
vec LSystem2d::get_as(){return avec;}
vec LSystem2d::get_bs(){return bvec;}

// Calculates rdot according to the (r,t) ODE
double LSystem2d::rdot(double r, vec * nvec)
{
	// nvec should be a unit vector, but normalize anyway
	double normn=norm(*nvec,2);	
	if ( normn < 1E-15)
		throw 5;
	vec nnvec = (*nvec) / normn;

	// Calculate (r,t) ODE
	return dot(bvec, nnvec) - r *(sum(avec) - avec(1)*nnvec(1)*nnvec(1) - 
		avec(2)*nnvec(2)*nnvec(2) - avec(0)*nnvec(0)*nnvec(0));
}

// This function does the heavy work for the thread: calculating fM(r), fm(r), nMVec(r) and nmvec(r) 
vec LSystem2d::threads(int mesh, string file)
{
	// Create .dat files to record fM(r)/fm(r), nvec(r)'s and the errors
	string filenamef = file + "_f.dat";
	string filenamen = file + "_n.dat";
	string filenamee = file + "_e.dat";
	fstream fdat, ndat, edat;

	fdat.open(filenamef.c_str(), fstream::out); 
	ndat.open(filenamen.c_str(), fstream::out); 
	edat.open(filenamee.c_str(), fstream::out);

	// Initialization of some parameters
	double stepsize = 1./(double)mesh;	
	double normb = norm(bvec,2);
	double horizon = 0; 
	vec nhorizon;

	// There are two threads: one each for nvec(0) = +/- bvec
	for (int m = 1; m >= -1; m -= 2){

		// We can solve exactly at r=0
		double rd = m*normb;
		ndat << 0 << " " << 0 << " " << 0 << endl;
		fdat << 0. << " " << rd << endl;
		edat << 0. << " " << 0. << endl;

		// Initial condition of (nvec,r)-ODE
		double r = 0; 
		vec nvec = m*bvec / normb; 

		double lambda, poly;
		for (int step = 1; step <= mesh; step++){
			
			// Apply Runge-Kutta
			new_RK(&nvec, r, stepsize); 
			r += stepsize;

			// Check to see if horizon if crossed
			rd = rdot(r, &nvec);
			if (rd > 0. && m==1) { 
				horizon = r; 
				nhorizon = nvec;
			}

			// Compute thread vector
			vec rnvec = r*rotation.t()*nvec;

			// Print updated variables for _f and _n files
			ndat 	<< rnvec(0) << " " << rnvec(1) << " " 
				<< rnvec(2) << " " << endl;
			fdat << r << " " << rd << endl;

			// Compute the error and print to _e file
			edat << r << " " << norm(cross(nvec, bvec+2*r*diagmat(avec)*nvec),2) << endl;
		}

		// Need blank lines to indicate finished loop
		ndat << endl; fdat << endl; edat << endl;
	}
	
	// Compute the horizon vector
	vec rnhor = horizon*rotation.t()*nhorizon;
	
	// Here follows the special cases where new threads may arise

	// Special case 1: a1 > a2 and b1 = 0
	if (bvec(0) == 0 && avec(0) != avec(1)){

		// Compute source of new thread
		double c1 = bvec(1)/2/(avec(0)-avec(1));
		double c2 = bvec(2)/2/(avec(0)-avec(2));
		double rstar = sqrt(pow(c1, 2) + pow(c2, 2));

		// Only matters if rstar < 1
		if (rstar < 1){
			
			// Compute one branch of new thread (backward to source)	
			for (int k = mesh; k >= ceil(rstar*(double)mesh); k-- ) {
				
				// Compute variables
				double r = (double)k/(double)mesh; 
				vec nvec("0 0 0");
				nvec(0) = sqrt(1-pow(rstar/r,2)); 
				nvec(1) = c1/r; 
				nvec(2) = c2/r;
				vec rnvec = r*rotation.t()*nvec;

				// Print to _n and _f
				ndat << rnvec(0) << " " << rnvec(1) << " " 
					<< rnvec(2) << endl;
				fdat << r << " " << rdot(r, &nvec) << endl;
			}

			// Compute second branch forward from source (so gnuplot
			//	sees one long thread)
			for (int k = ceil(rstar*(double)mesh); k <= mesh; k++ ) {
				double r = (double)k/(double)mesh; 
				vec nvec("0 0 0");
				nvec(0) = -sqrt(1-pow(rstar/r,2)); 
				nvec(1) = c1/r; 
				nvec(2) = c2/r;
				vec rnvec = r*rotation.t()*nvec;
				ndat << rnvec(0) << " " << rnvec(1) << " " 
					<< rnvec(2) << endl;
				fdat << r << " " << rdot(r, &nvec) << endl;
			}
			
			// Blank lines to indicate thread is finished
			fdat << endl; 
			ndat << endl;
		}
	}

	// Add dummy line to _f file since gnuplot will complain if not there
	else 	
		fdat << 0 << " " << 0 << endl << endl;
		
	// Special case 2: a1 > a2 > a3 and b2 = 0
	if (bvec(1) == 0 && avec(1) != avec(0) && avec(1) != avec(2)){

		// Compute source of new thread
		double c1 = bvec(0)/2/(avec(1)-avec(0)); 
		double c2 = bvec(2)/2/(avec(1)-avec(2));
		double rstar = sqrt(pow(c1, 2) + pow(c2, 2));

		// Only matters if rstar <1
		if (rstar < 1){
			
			// Compute one branch of new thread (backward to source)	
			for (int k = mesh; k >= ceil(rstar*(double)mesh); k-- ) 
			{
				// Compute variables
				double r = (double)k/(double)mesh; 
				vec nvec("0 0 0");
				nvec(1) = sqrt(1-pow(rstar/r,2)); 
				nvec(0) = c1/r; 
				nvec(2) = c2/r;
				vec rnvec = r*rotation.t()*nvec;
				
				// Print to _n and _f
				ndat << rnvec(0) << " " << rnvec(1) << " " 
					<< rnvec(2) << endl;
				fdat << r << " " << rdot(r, &nvec) << endl;
			}
		
			// Compute second branch forward from source (so gnuplot
			//	sees one long thread)
			for (int k = ceil(rstar*(double)mesh); k <= mesh; k++ ) {

				double r = (double)k/(double)mesh; 
				vec nvec("0 0 0");
				nvec(1) = -sqrt(1-pow(rstar/r,2)); 
				nvec(0) = c1/r; 
				nvec(2) = c2/r;
				vec rnvec = r*rotation.t()*nvec;

				ndat << rnvec(0) << " " << rnvec(1) << " " 
					<< rnvec(2) << endl;
				fdat << r << " " << rdot(r, &nvec) << endl;
			}

			// Blank lines to indicate thread is finished
			fdat << endl; 
			ndat << endl;
		}
	}
	// Add dummy line to _f file since gnuplot will complain if not there
	else 	
		fdat << 0 << " " << 0 << endl << endl;


	// Special case 3: a2 > a3 and b3 = 0
	if (bvec(2) == 0 && avec(2) != avec(0) && avec(2) != avec(1)){

		// Compute source of new thread
		double c1 = bvec(0)/2/(avec(2)-avec(0)); 
		double c2 = bvec(1)/2/(avec(2)-avec(1));
		double rstar = sqrt(pow(c1, 2) + pow(c2, 2));

		// Only matters if rstar < 1
		if (rstar < 1){
		
			// Compute one branch of new thread (backward to source)	
			for (int k = mesh; k >= ceil(rstar*(double)mesh); k-- ) {

				// Compute variables
				double r = (double)k/(double)mesh; 
				vec nvec("0 0 0");
				nvec(2) = sqrt(1-pow(rstar/r,2)); 
				nvec(0)=c1/r; 
				nvec(1) = c2/r;
				vec rnvec = r*rotation.t()*nvec;
								
				// Print new values to _f and _n files
				ndat << rnvec(0) << " " << rnvec(1) << " " 
					<< rnvec(2) << endl;
				fdat << r << " " << rdot(r, &nvec) << endl;
			}

			// Compute second branch forward from source (so gnuplot
			//	sees one long thread)
			for (int k = ceil(rstar*(double)mesh); k <= mesh; k++ ) {
				double r = (double)k/(double)mesh; 
				vec nvec("0 0 0");
				nvec(2) = -sqrt(1-pow(rstar/r,2)); 
				nvec(0) = c1/r; 
				nvec(1) = c2/r;
				vec rnvec = r*rotation.t()*nvec;

				// Print 
				ndat << rnvec(0) << " " << rnvec(1) << " " 
					<< rnvec(2) << endl;
				fdat << r << " " << rdot(r, &nvec) << endl;
			}

			fdat << endl; 
			ndat << endl;
		}
	}

	// Special case 4: a1 = a2 and b1 = b2 = 0
	if (avec(0) == avec(1) && bvec(0) == 0 && bvec(1) == 0){

		// Compute source of new threads
		double c = bvec(2)/2/(avec(0)-avec(2)); 
		
		// Only matters if c < 1
		if (abs(c) < 1){
			
			// Loop over spokes
			for (int m = 0; m < SPOKES; m++){

				// Calculate spoke ending
				double c12 = sqrt(1-c*c);
				double c1 = c12*cos(2*PI*m/SPOKES);
				double c2 = c12*sin(2*PI*m/SPOKES);	

				
				vec cvec("0 0 0");
				cvec(2) = c;
				
				// Loop over points in spoke
				for (int k = 0; k <= mesh; k++ ) {

					// Calculate point
					double medc1 = c1*(1-abs(2*k/mesh-1));
					double medc2 = c2*(1-abs(2*k/mesh-1));
					cvec(0) = medc1;
					cvec(1) = medc2;
					vec rnvec = rotation.t()*cvec;
				
					// Print new values to _n file
					ndat << rnvec(0) << " " << rnvec(1) << " " << rnvec(2) << endl;
				}		
			}
			
			// Calculate rdot and print
			for (int k = ceil(c*(double)mesh); k <= mesh; k++) {
				double r = k/mesh;
				double fr = -(avec(0) + avec(2))*r + bvec(2)*bvec(2)/4/r/(avec(0) - avec(2));
				fdat << r << " " << fr << endl;
			}
			fdat << endl; 
			ndat << endl;
		}
	}

	// Special case 5: a2 = a3 and b2 = b3 = 0
	if (avec(1) == avec(2) && bvec(1) == 0 && bvec(2) == 0){

		// Compute source of new threads
		double c = bvec(0)/2/(avec(2)-avec(0)); 
		
		// Only matters if c < 1
		if (abs(c) < 1){
			
			// Loop over spokes
			for (int m = 0; m < SPOKES; m++){

				// Calculate spoke ending
				double c12 = sqrt(1-c*c);
				double c1 = c12*cos(2*PI*m/SPOKES);
				double c2 = c12*sin(2*PI*m/SPOKES);	

				
				vec cvec("0 0 0");
				cvec(0) = c;
				
				// Loop over points in spoke
				for (int k = 0; k <= mesh; k++ ) {

					// Calculate point
					double medc1 = c1*(1-abs(2*k/mesh-1));
					double medc2 = c2*(1-abs(2*k/mesh-1));
					cvec(1) = medc1;
					cvec(2) = medc2;
					vec rnvec = rotation.t()*cvec;
				
					// Print new values to _n file
					ndat << rnvec(0) << " " << rnvec(1) << " " << rnvec(2) << endl;
				}		
			}
			
			// Calculate rdot and print
			for (int k = ceil(c*(double)mesh); k <= mesh; k++) {
				double r = k/mesh;
				double fr = -(avec(0) + avec(2))*r + bvec(0)*bvec(2)/4/r/(avec(2) - avec(0));
				fdat << r << " " << fr << endl;
			}
			fdat << endl; 
			ndat << endl;
		}
	}
	
	// Add dummy line to _f file since gnuplot will complain if not there
	else 	
		fdat << 0 << " " << 0 << endl << endl;
		
	// Close data files and return the horizon value (non-special version)
	fdat.close(); ndat.close();edat.close();
	return rnhor;
}

// For a given (nvec, r) on S1 cross [0,1], computes d/dt(nvec) according to
//	the (nvec,r) ODE
vec LSystem2d::ndot(double r, vec* nvec){
	
	// Enforce normalization
	vec nv = *nvec / norm(*nvec,2); 

	// Define projection matrix perp to nvec 
	mat Porth = eye(3,3) - nv*nv.t();

	// Helper matrices
	mat AS = diagmat(avec); 
	double Lambda = dot(bvec,nv)+2*r*dot(nv, AS*nv);

	// Compute in the case that 2rAS-Lambda is not invertible
	if ( 	(avec(0) == 0 && nv(0) == 0) || 
		(avec(1) == 0 && nv(1) == 0) || 
		(avec(2) == 0 && nv(2) == 0)  )
		return 2*Porth*AS*nv/norm(Porth*(2*r*AS-Lambda*eye(3,3)),2);
	
	// Compute in the invertible case
	else{
		mat rASLi = inv(2*r*AS-Lambda*eye(3,3));
      		double k = -2*dot(nv,rASLi*AS*nv)/dot(nv,rASLi*nv);
		return rASLi*(-k*eye(3,3)-2*AS)*nv;
	}
}

// This function computes the chimney
void LSystem2d::chimney(int mesh, string file, vec rnhor) {

	// Initialization
	double stepsize = 1. / (double)mesh;
	vec wall = zeros<vec>(3);

	// Open new data files for chimney and errors
	string filenamew = file + "_w.dat";
	string filenamewe = file + "_we.dat";
	fstream wdat, wedat;
	wdat.open(filenamew.c_str(), fstream::out);
	wedat.open(filenamewe.c_str(), fstream::out);

	// Sweep over the lines in the chimney
	for (int m = 0; m < CHIMNEY_LINES; m++) {

		// Compute start of the chimney line
		double r = 0;
		double angle = m * 2 * PI / CHIMNEY_LINES;
		vec unit1 = zeros<vec>(3);
		unit1(1) = bvec(2);
		unit1(2) = -bvec(1);
		if (norm(unit1) <= 1e-10) {
			unit1(1) = 1;
		}
		unit1 = normalise(unit1);
		vec unit2 = cross(bvec, unit1);
		unit2 = normalise(unit2);
		wall = cos(angle)*unit1 + sin(angle)*unit2;
		
		// Record initial info
		vec rwvec = r*rotation.t()*wall;
		wdat << rwvec(0) << " " << rwvec(1) << " "
			<< rwvec(2) << endl;
		wedat << 0 << " " << 0 << endl;

		// Initialize the exit parameter. Can't start at zero, else chimney
		// will terminate prematurely
		double exit_param = 100.;

		do {
			// Compute next point in chimney line
			vec oldwall = wall;
			new_wall_RK(&wall, r, stepsize);
			wall = normalise(wall);
			rwvec = r*rotation.t()*wall;

			// Print updated variables for _w file
			wdat << rwvec(0) << " " << rwvec(1) << " "
				<< rwvec(2) << " " << endl;

			// Increment and compute error
			r += stepsize;
			double werr = abs(rdot(r, &wall));
			wedat << r << " " << werr << endl;

			// Compute exit parameter
			vec grad = bvec + 2 * r*(avec % wall);
			double Lambda = dot(wall, grad);			
			exit_param = abs(dot(grad, grad) / Lambda / Lambda - 1);

			// Unfortunately algorithm is finicky near the exit. 
			// Threshold error can't be too small, but if too big will look
			// articially flat at the end. Also must insert final value
			// by hand.
		} while (exit_param > CHIMNEY_TOL);

		// In one special case, the final values are not necessarily the 
		// typical horizon.
		double rn2t = bvec(1) / 2 / (avec(0) - avec(1));
		double rn3t = bvec(2) / 2 / (avec(0) - avec(2));
		double thresh = sqrt(rn2t*rn2t + rn3t*rn3t);
		if (bvec(0) == 0 && avec(0) > avec(1) && thresh < norm(rnhor)) {
			
			// Calculate alternate horizons
			vec alt_rnhor = zeros<vec>(3);
			double rsq = bvec(1) * bvec(1) / 4 / (avec(0) - avec(1)) / (avec(1) + avec(2));
			rsq += bvec(2) * bvec(2) / 4 / (avec(0) - avec(2)) / (avec(1) + avec(2));

			
			// Exactly which horizon to insert depends on which line
			if (m > 0 && m < CHIMNEY_LINES/2 && rsq <= 1) {
				alt_rnhor(1) = bvec(1) / 2 / (avec(0) - avec(1));
				alt_rnhor(2) = bvec(2) / 2 / (avec(0) - avec(2));
				alt_rnhor(0) = - sqrt(rsq - alt_rnhor(1) * alt_rnhor(1) - alt_rnhor(2) * alt_rnhor(2));
			}
			else if (m > CHIMNEY_LINES/2 && rsq <= 1) {
				alt_rnhor(1) = bvec(1) / 2 / (avec(0) - avec(1));
				alt_rnhor(2) = bvec(2) / 2 / (avec(0) - avec(2));
				alt_rnhor(0) = sqrt(rsq - alt_rnhor(1) * alt_rnhor(1) - alt_rnhor(2) * alt_rnhor(2));
			}
			else alt_rnhor = rnhor;
			alt_rnhor = rotation.t()*alt_rnhor;

			// Insert final value ten times since gnuplot will skip values, and we don't want
			// the horizons skipped.
			if (norm(alt_rnhor - rwvec) < norm(rnhor - rwvec)){
				rnhor = alt_rnhor;
			}	
			if (norm(rnhor - rwvec) < 0.1){
				for (int p = 1; p <= 10; p++) {
					wdat << alt_rnhor(0) << " " << alt_rnhor(1) << " " << alt_rnhor(2) << endl;
					wedat << norm(alt_rnhor) << " " << 0 << endl;
				}
			}
		}

		// Outside of the special case, the final value is often the horizon calculated in the threads 
		// function. If the last point in the chimney line is not sufficiently close to this horizon,
		// the final value is omitted.
		else if (norm(rnhor - rwvec) < 0.1){
			for (int p = 1; p <= 10; p++) {
				wdat << rnhor(0) << " " << rnhor(1) << " " << rnhor(2) << endl;
				wedat << norm(rnhor) << " " << 0 << endl;
			}
		}

		// Empty lines to indicate finish of a line
		wdat << endl;
		wedat << endl;
	}

	wdat.close(); wedat.close();
}

// This function computes the correction to stay on the chimney
vec LSystem2d::chimneyCorr(double r, vec* wall) {

	// Compute auxiliary quantities
	vec grad = bvec + 2 * r*(avec % (*wall));
	double Lambda = dot(*wall, grad);
	double numer = sum(avec - avec % *wall % *wall);
	double denom = dot(grad, grad) - Lambda*Lambda;

	// If denom is zero, we reached the chimneytop
	if (denom == 0.)
		return zeros<vec>(3);

	// Calculate correction
	else
		return (grad + (*wall)*Lambda)*numer / denom;
}


// Creates the output files for gnuplot
void LSystem2d::gnufilescreate(string filename)
{
	fstream plotfile;

	// Create "***_f_plot.gnu" which plots fM(r) and fm(r) 
	string filenamef = filename + "_f_plot.gnu";	
	plotfile.open(filenamef.c_str(), fstream::out);
	plotfile << "reset" << endl << endl
		<< "set terminal windows size 700,524 enhanced font 'Verdana,10'" 
		<< endl 
		<< "set border lw 1.5; unset key; set view 53,16" << endl
		<< "set xrange [-0:1]; set xlabel \"r\"; set ylabel \"dr/dt\"; f(x) = 0"
		<< endl << endl

		<< "plot f(x) lc rgb \"black\" lw 0.5, \"" << filename 
		<< "_f.dat\" every :::0::0 lt 6 lw 2 lc rgb \"black\" with lines, \""
		<< filename 
		<< "_f.dat\" every :::1::1 lt 6 lw 2 lc rgb \"black\" with lines";

	// For some input, there are multiple branches, so 
		plotfile << ", \"" << filename << "_f.dat\" every 10:::2::2 "
			"lt 0 lw 2 lc rgb \"black\" with lines, \"" << filename 
		<< "_f.dat\" every 10:::3::3 lt 0 lw 2 lc rgb \"black\" with lines, \""
		<< filename << "_f.dat\" every 10:::4::4 lt 0 lw 2 lc rgb \"black\""
		<< "with lines" << endl;
	
	plotfile << endl;
	plotfile.close();

	// Make "***_e_plot.gnu" which plots the error in the thread
	string filenamee = filename + "_e_plot.gnu";
	plotfile.open(filenamee.c_str(), fstream::out);
	plotfile << "reset" << endl << endl
		<< "set terminal windows size 700,524 enhanced font 'Verdana,10'" 
		<< endl
		<< "set border lw 1.5; unset key; set view 53,16" << endl
		<< "set xrange [-0:1]; set xlabel \"r\";" 
			"set ylabel \"||D_n f(n,r)||\"; f(x) = 0" << endl << endl
		<< "plot \"" << filename 
		<< "_e.dat\" every :::0::0 lt 6 lw 2 lc rgb \"black\" with lines, \""
		<< filename 
		<< "_e.dat\" every :::1::1 lt 6 lw 2 lc rgb \"black\" with lines" 
		<< endl;
	plotfile.close();

	// Make "***_we_plot.gnu" which plots the error in the chimney
	string filenamewe = filename + "_we_plot.gnu";
	plotfile.open(filenamewe.c_str(), fstream::out);
	plotfile << "reset" << endl << endl
		<< "set terminal windows size 700,524 enhanced font 'Verdana,10' "
		<< endl
		<< "set border lw 1.5; unset key; set view 53,16" << endl
		<< "set xrange [-0:1]; set xlabel \"r\";"
		"set ylabel \"||f(n,r)||\"; f(x) = 0" << endl << endl
		<< "plot \"" << filename
		<< "_we.dat\" every :::0::0 lt 6 lw 2 lc rgb \"black\" with lines";
	for (int m = 1; m < CHIMNEY_LINES; m++) {
		plotfile << ", \"" << filename << "_we.dat\" every :::" << m << "::" << m
			<< " lt 6 lw 2 lc rgb \"black\" with lines";
	}
	plotfile.close();

	// Make "***_n_plot.gnu" which plots the thread and chimney
	string filenamen = filename + "_n_plot.gnu";
	plotfile.open(filenamen.c_str(), fstream::out);
	plotfile << "reset" << endl << endl
		<< "set termoption dashed" << endl
		<< "set terminal windows size 700,524 enhanced font 'Verdana,10'"
		<< endl
		<< "set border lw 1.5" << endl
		<< "unset key; unset tics; unset border" << endl
		<< "set lmargin screen 0" << endl
		<< "set rmargin screen 1" << endl
		<< "set tmargin screen 1" << endl
		<< "set bmargin screen 0" << endl
		<< "set size ratio -1" << endl
		<< "set view 53,16" << endl
		<< "set parametric" << endl
		<< "set isosamples 30" << endl
		<< "set hidden3d" << endl
		<< "set xrange [-1.2:1.2]" << endl
		<< "set yrange [-1.2:1.2]" << endl
		<< "set zrange [-1.2:1.2]" << endl
		<< "set urange [0:3.0/2*pi]" << endl
		<< "set vrange [-pi/2:pi/2]" << endl << endl
		<< "r = 1.0" << endl
		<< "fx(v,u) = r*cos(v)*cos(u)" << endl
		<< "fy(v,u) = r*cos(v)*sin(u)" << endl
		<< "fz(v)   = r*sin(v)" << endl << endl
		<< "splot fx(v,u),fy(v,u),fz(v) lc rgb \"grey\" lw 0.5, \""
		<< filename << "_n.dat\" every 10:::0::0 lt 6 lw 2 "
		"lc rgb \"black\" with lines, \""
		<< filename << "_n.dat\" every 10:::1::1 lt 0 lw 2 lc rgb"
		" \"black\" with lines, \""
		<< filename << "_n.dat\" every 100:::2::2 lt 10 lw 1 lc rgb"
		" \"black\" with lines";

	for (int m = 0; m < CHIMNEY_LINES; m++) {
		plotfile << ", \"" << filename << "_w.dat\" every 6:::" << m << "::" << m << " lt 7 lw 1 lc rgb \"black\" with lines"; 
	}
		plotfile << endl;
	plotfile.close();
}


