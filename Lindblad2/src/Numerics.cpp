/* 
 * Lindblad2
 *
 * Copyright (c) 2014 - 2015 Patrick Rooney (darraghrooney@gmail.com)
 *
 * This file is part of Lindblad2.
 *
 * Lindblad2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Lindblad2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Lindblad2.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "LSystem2d.h"

using namespace std;
using namespace arma;

// This function computes one step of the (nvec,r) ODE
// 	using Euler's method (on manifolds)
void LSystem2d::newn_Euler(vec* nvec, double r, double dr){
	
	// Compute d/dt(nvec) at r, and decompose
	vec nd = ndot(nvec, r);
	double theta = norm(nd,2); 
	vec ndn = nd/theta;
	
	// Apply Euler method for this manifold
	if (theta != 0) 
		*nvec = (*nvec)*cos(dr*theta) + ndn*sin(dr*theta);
}

// This function computes one step of the numerical evolution on So(3) x [0,1]
// 	using Heun's method (on manifolds)
void LSystem2d::newn_Heun(vec* nvec, double r, double dr){
	
	// Compute k1
	vec zerovec = vec("0 0 0");
	vec k1 = newk(nvec, &zerovec, r, dr);
	
	if (norm(k1,2) != 0){
		
		// Compute k2 and ktot		
		vec k2 = newk(nvec, &k1, r+dr, dr);
		vec ktot = (k1+k2)/2; 

		// Apply Heun's method
		double theta = norm(ktot,2);
		vec ktotn = ktot/theta;
		*nvec = (*nvec)*cos(theta) + cross(ktotn,*nvec)*sin(theta);}
}


// This function computes one step of the numerical evolution on So(3) x [0,1]
// 	using the Runge-Kutta method (on manifolds)
void LSystem2d::newn_RK(vec* nvec, double r, double dr){

	// Compute k1
	vec zerovec = vec("0 0 0");
	vec k1 = newk(nvec, &zerovec, r, dr);

	if (norm(k1,2) != 0){

		// Compute the other k's
		vec k1i = k1/2.;
		vec k2 = newk(nvec, &k1i, r+dr/2., dr);
		vec k2i = k2/2. - cross(k1,k2)/8.;
		vec k3 = newk(nvec, &k2i, r+dr/2, dr);
		vec k4 = newk(nvec, &k3, r+dr, dr);
		vec ktot = (k1 + 2*k2 + 2*k3 + k4 - cross(k1,k4)/2.)/6.; 
		
		// Apply Runge-Kutta
		double theta = norm(ktot,2);
		vec ktotn = ktot/theta;
		*nvec = (*nvec)*cos(theta) + cross(ktotn,*nvec)*sin(theta);}
}

// This is a helper function for the above numerical methods
vec LSystem2d::newk(vec* nvec, vec* oldk, double r, double dr){

	// Case for non-zero oldk
	double theta = norm(*oldk,2);
	if (theta != 0) {

		// Compute helper quantities
		vec oldkn = *oldk/theta; 
		vec ndn = -cross(*nvec,oldkn);
		vec nt = cos(theta)*(*nvec)+sin(theta)*ndn;
		vec ki = cross(nt,ndot(&nt, r));

		// Compute newk
		vec newk = dr*dot(ki,oldkn)*oldkn + dr*(dot(ki,*nvec)*cos(theta)
			+ dot(ki,ndn)*sin(theta))*(*nvec)+dr*(-dot(ki,*nvec)*sin(theta)
			+ dot(ki,ndn)*cos(theta))*ndn;
		return newk;
	}

	// Case for zero oldk
	else
		return dr*cross(*nvec, ndot(nvec,r));
}

