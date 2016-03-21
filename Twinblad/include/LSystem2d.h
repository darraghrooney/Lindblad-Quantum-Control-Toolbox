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


#ifndef LSYSTEM2D_H_INCLUDED
#define LSYSTEM2D_H_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <ctime>
#include <functional>

#define PI 3.141592653589793238462643383

using namespace std;
using namespace arma;

class LSystem2d
{

private:
	// Parameters specifying the GKS matrix 
	vec avec, bvec;								
	mat rotation;								

public:
	// Default constructor has avec("1 1 0") and bvec("0 0 2").
	LSystem2d(void);							
	// This constructor takes avec and bvec. Rotation is the identity.
	LSystem2d(vec *, vec *);
	// This construction takes avec, bvec and a rotation matrix.
	LSystem2d(vec *, vec *, mat *);			
	// This constructs a random system if 'r' is given.
	LSystem2d(char);				
	
	// Accesses the GKS parameters
	vec get_as();			
	vec get_bs();

	// Computes dr/dt = f(r,nvec) using avec and bvec
	double rdot(double, vec*);		

	// Computes the time derivative of nvec given nvec and r
	vec ndot(double, vec*);	

	// This is the main function that computes the threads
	vec threads(int, string);

	// This is the function that computes the chimney
	void chimney(int, string, vec);
	vec chimneyCorr(double, vec*);

	// Constructs the .gnu files that plot the f, n and e data
	void gnufilescreate(string);

	// These are helper functions for the numerics
	void new_Euler(vec*, double, double);
	void new_Heun(vec*, double, double);
	void new_RK(vec*, double, double);
	vec newk(vec*, vec*, double, double);
	
	// Numeric versions intended for the chimney
	void new_wall_RK(vec*, double, double);
	vec new_wall_k(vec*, vec*, double, double);

};

#endif // LSYSTEM2D_H_INCLUDED
