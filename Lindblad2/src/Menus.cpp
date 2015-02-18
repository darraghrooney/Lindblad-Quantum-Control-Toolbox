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

#include <iostream>
#include <vector>
#include <complex>
#include <armadillo>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "Menus.h"

using namespace std;
using namespace arma;

char mainMenu(void){
	
	// Secondary menu is chosen based on a character input from console
	char choice;

	cout << endl << "/*** Lindblad2 Copyright (C) 2015 Patrick Rooney ***/" << endl << endl;
    	cout << "This program comes with ABSOLUTELY NO WARRANTY. This is free software, "
		"and you are welcome to redistribute it under certain conditions. "
		"For details, consult the GNU Personal License "
		"that came with this software. " << endl << endl;

    	cout << "Welcome to Lindblad2. " << endl;
	cout << "This is a program for optimizing quantum control of "
			"n=2 Lindblad systems." << endl << endl;
	cout << "Press 'r' for random GKS matrix." << endl;
	cout << "Press 'f' to enter data as a file." << endl;
	cout << "Press 's' to enter data for the real-diagonalized Pauli GKS matrix." << endl;
	cout << "Press 'p' to enter data for the undiagonalized Pauli GKS matrix." << endl;
	cout << "Press 'j' to enter data for the +/-/z GKS matrix." << endl << endl;
	cin >> choice; cout << endl;
	
	return choice;
}

// Menu for randomized input (input here being a_j, b_j, j=1,2,3)
void randomMenu(void){
	int reps; 
	string file;

	// Construct random GKS system
	LSystem2d randsys('r');

	// Print a_j, b_j to console	
	cout << "Parameters of random GKS matrix are: " << endl << endl;
	cout << "a1 = " << randsys.get_as().at(0) << endl;
	cout << "a2 = " << randsys.get_as().at(1) << endl;
	cout << "a3 = " << randsys.get_as().at(2) << endl;
	cout << "b1 = " << randsys.get_bs().at(0) << endl;
	cout << "b2 = " << randsys.get_bs().at(1) << endl;
	cout << "b3 = " << randsys.get_bs().at(2) << endl << endl;
	
	// Prompt for #reps and output filenames 
	cout << "Enter no. reps (1000 recommended): "; cin >> reps;
	cout << "Enter filename (without extension): "; cin >> file;

	// Create parameter file
	fstream paramfile;

	string pfile = file + "_param.txt";	
	paramfile.open(pfile.c_str(), fstream::out);
	paramfile << "a1 = " << randsys.get_as().at(0) << endl;
	paramfile << "a2 = " << randsys.get_as().at(1) << endl;
	paramfile << "a3 = " << randsys.get_as().at(2) << endl;
	paramfile << "b1 = " << randsys.get_bs().at(0) << endl;
	paramfile << "b2 = " << randsys.get_bs().at(1) << endl;
	paramfile << "b3 = " << randsys.get_bs().at(2) << endl << endl;
	paramfile.close();
	
	// Do the calculations and create .dat output files

	cout << "Running algorithm..." << endl;
	randsys.threads(reps, file);

	// Create the .gnu files
	randsys.gnufilescreate(file) ; 
}

// Menu for input from file
void fileMenu(void){
	ifstream input; 
	string infile, outfile;
	int reps; 
	
	// Instructions for input file format (need six real numbers)
    	cout << "The input file should contain six real numbers." << endl
		<< "They will be assigned as a1, a2, a3, b1, b2, b3, in that order." << endl
		<< "All data after the sixth entry will be disregarded. " << endl
		<< "If there are fewer than six, the remaining variables"
		" will be set to zero." << endl << endl;
    	cout << "Enter filename, including extension: " ;
    	cin >> infile;
	
	// Re-prompt if input file doesn't exist
	while (!ifstream(infile.c_str())){
		cout << "File doesn't exist." << endl;
	    	cout << "Enter filename, including extension: " ;
	    	cin >> infile;
	}
	
    	vec as(3); 
	vec bs(3);
	
	// Read input and print to console
	input.open(infile.c_str(), fstream::in);
	input >> as(0) >> as(1) >> as(2) >> bs(0) >> bs(1) >> bs(2);
	cout << "a1 = " << as(0) << endl;
	cout << "a2 = " << as(1) << endl;
	cout << "a3 = " << as(2) << endl;
	cout << "b1 = " << bs(0) << endl;
	cout << "b2 = " << bs(1) << endl;
	cout << "b3 = " << bs(2) << endl << endl;
	input.close();

	// Construct the Lindblad system
	LSystem2d insys(&as, &bs);;
	
	// Ask for #reps and output filename
	cout << "Enter no. reps (1000 recommended): "; cin >> reps;
	cout << "Enter filename (without extension): "; cin >> outfile;
	cout << "Running algorithm..." << endl;

	// Do the calculations and create .dat output files
	insys.threads(reps, outfile);
	// Create the .gnu files
	insys.gnufilescreate(outfile); 

}

// Menu for console input, where input is GKS matrix in Pauli basis
void PauliMenu(void){
	int reps; 
	string file; 
	vec diag(3), reoff(3), imoff(3);

	// Ask for GKS elements
	cout << "Enter axx: "; cin >> diag(0);
	cout << "Enter ayy: "; cin >> diag(1);
	cout << "Enter azz: "; cin >> diag(2);
	cout << "Enter Re ayz: "; cin >> reoff(0);
	cout << "Enter Re azx: "; cin >> reoff(1);
	cout << "Enter Re axy: "; cin >> reoff(2);
	cout << "Enter Im ayz: "; cin >> imoff(0);
	cout << "Enter Im azx: "; cin >> imoff(1);
	cout << "Enter Im axy: "; cin >> imoff(2);

	// Construct real component of the GKS matrix from input
	mat gkssym(3,3); 
	gkssym.diag() = diag;
	gkssym(0,1) = reoff(0); 
	gkssym(1,0) = reoff(0);
	gkssym(0,2) = reoff(1); 
	gkssym(2,0) = reoff(1);
	gkssym(1,2) = reoff(2); 
	gkssym(2,1) = reoff(2);
	
	// Calculate rotation matrix to diagonalize gkssym
	vec as(3), bs(3); 
	mat rotation(3,3);
        eig_sym(as,rotation,gkssym);

	// Determine avec and bvec
	as.swap_rows(0,2); 
	rotation.swap_cols(0,2);
	bs = 2*rotation.t()*imoff;

	// Print a_j and b_j to console
	cout << endl << "a1 = " << as(0) << endl;
	cout << "a2 = " << as(1) << endl;
	cout << "a3 = " << as(2) << endl;
	cout << "b1 = " << bs(0) << endl;
	cout << "b2 = " << bs(1) << endl;
	cout << "b3 = " << bs(2) << endl << endl;
	cout << "Rotation matrix: " << endl << rotation << endl << endl;

	// Construct the Lindblad system
	LSystem2d paulisys(&as, &bs, &rotation);

	// Ask for #reps and output filename
	cout << "Enter no. reps (1000 recommended): "; cin >> reps;
	cout << "Enter filename (without extension): "; cin >> file;

	// Create parameter file
	fstream paramfile;

	string pfile = file + "_param.txt";	
	paramfile.open(pfile.c_str(), fstream::out);
	paramfile << "a1 = " << as(0) << endl;
	paramfile << "a2 = " << as(1) << endl;
	paramfile << "a3 = " << as(2) << endl;
	paramfile << "b1 = " << bs(0) << endl;
	paramfile << "b2 = " << bs(1) << endl;
	paramfile << "b3 = " << bs(2) << endl << endl;
	paramfile << "Rotation matrix: " << endl << rotation << endl;
	paramfile.close();

	// Do the calculations and create .dat output files
	cout << "Running algorithm..." << endl;
	paulisys.threads(reps, file);

	// Create the .gnu files
	paulisys.gnufilescreate(file);
}

// Menu for console input, where input is the avec and bvec
void realdiagMenu(void){
	int reps; 
	string file; 
	vec as(3), bs(3);

	// Ask for input
	cout << "Enter a1: "; cin >> as(0);
	cout << "Enter a2: "; cin >> as(1);
	cout << "Enter a3: "; cin >> as(2);
	cout << "Enter b1: "; cin >> bs(0);
	cout << "Enter b2: "; cin >> bs(1);
	cout << "Enter b3: "; cin >> bs(2);

	// Construct Lindblad system
	LSystem2d specsys(&as, &bs);

	// Ask for #reps and output filename
	cout << "Enter no. reps (1000 recommended): "; cin >> reps;
	cout << "Enter filename (without extension): "; cin >> file;

	// Create parameter file
	fstream paramfile;

	string pfile = file + "_param.txt";	
	paramfile.open(pfile.c_str(), fstream::out);
	paramfile << "a1 = " << as(0) << endl;
	paramfile << "a2 = " << as(1) << endl;
	paramfile << "a3 = " << as(2) << endl;
	paramfile << "b1 = " << bs(0) << endl;
	paramfile << "b2 = " << bs(1) << endl;
	paramfile << "b3 = " << bs(2) << endl << endl;
	paramfile.close();

	// Do calculations and create data output files
	cout << "Running algorithm..." << endl;
	specsys.threads(reps, file);
	// Create files for gnuplot
	specsys.gnufilescreate(file);
}

// Similar to PauliMenu, but basis is +/-/z instead of x/y/z
void pmzMenu(void){
	int reps; 
	string file; 
	vec gks(9), diag(3), reoff(3), imoff(3);

	// Ask for GKS elements
	cout << "Enter a++: "; cin >> gks(0);
	cout << "Enter a--: "; cin >> gks(1);
	cout << "Enter azz: "; cin >> gks(2);
	cout << "Enter Re a+-: "; cin >> gks(3);
	cout << "Enter Re a+z: "; cin >> gks(4);
	cout << "Enter Re a-z: "; cin >> gks(5);
	cout << "Enter Im a+-: "; cin >> gks(6);
	cout << "Enter Im a+z: "; cin >> gks(7);
	cout << "Enter Im a-z: "; cin >> gks(8);

	// Calculate GKS matrix in Pauli basis
	diag(0) = gks(0)+gks(1)+gks(3)*2;
	diag(1) = gks(0)+gks(1)-gks(3)*2;
	diag(2) = gks(2);
	reoff(0)= gks(7)-gks(8);
	reoff(1)= gks(4) + gks(5);
	reoff(2)= -2 * gks(6);
	imoff(0)= - gks(4) + gks(5);
	imoff(1)= gks(7) + gks(8);
	imoff(2)= gks(0) - gks(1);

	// Construct symmetrized GKS matrix
	mat gkssym(3,3); 
	gkssym.diag() = diag;
	gkssym(0,1) = reoff(0); 
	gkssym(1,0) = reoff(0);
	gkssym(0,2) = reoff(1); 
	gkssym(2,0) = reoff(1);
	gkssym(1,2) = reoff(2); 
	gkssym(2,1) = reoff(2);

	// Diagonalize gkssym
    	vec as(3), bs(3); 
	mat rotation(3,3);
	eig_sym(as,rotation,gkssym);

	// Calculate a_j and b_j
	as.swap_rows(0,2); 
	rotation.swap_cols(0,2);
	bs = 2*rotation.t()*imoff;

	// Construct Lindblad system 
	LSystem2d paulisys(&as, &bs, &rotation);

	// Ask for #reps and output filename
	cout << "Enter no. reps: "; cin >> reps;
	cout << "Enter filename: "; cin >> file;

	// Create parameter file
	fstream paramfile;

	string pfile = file + "_param.txt";	
	paramfile.open(pfile.c_str(), fstream::out);
	paramfile << "a1 = " << as(0) << endl;
	paramfile << "a2 = " << as(1) << endl;
	paramfile << "a3 = " << as(2) << endl;
	paramfile << "b1 = " << bs(0) << endl;
	paramfile << "b2 = " << bs(1) << endl;
	paramfile << "b3 = " << bs(2) << endl << endl;
	paramfile << "Rotation matrix: " << endl << rotation << endl;
	paramfile.close();
	
	// Do calculations and create data output files
	cout << "Running algorithm..." << endl;
	paulisys.threads(reps, file);
	// Create files for gnuplot
	paulisys.gnufilescreate(file);
}
