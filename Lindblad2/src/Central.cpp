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
#include "Menus.h"
#include "LSystem2d.h"

#define PI 3.141592653589793238462643383

using namespace std;
using namespace arma;

int main()
{
	// Main function runs main menu then diverts to the secondary menu
	// 	based on choice of input format.	

	char choice = mainMenu();
	try{
		// Random input parameters
		if(choice == 'r')
			randomMenu();	

		// Input parameters from a file	
		else if (choice == 'f')
			fileMenu();
		
		// Console input wrt Pauli matrices
		else if (choice == 'p')
			PauliMenu();
		
		// Console input in avec, bvecv form 
		else if (choice == 's')
			realdiagMenu();
		
		// Console input in +/-/z basis form
		else if (choice == 'j')
			pmzMenu(); 
		
		else 	
			throw 0;

		cout << "Done. Files can now be processed with Gnuplot." << endl;
	}
	catch(int e){
		switch(e){
		
			// Thrown in mainMenu()
			case 0: cout << "You didn't enter an appropriate character."
					" Exiting program." << endl;
				break;
			
			// Cases 1-5 thrown in LSystem2d constructors (invalid input params)
			case 1:	cout << "The a_j's are not sorted or one is negative."	
					" Exiting program." << endl;
				break;
			case 2:	cout << "Positive semi-definiteness violated."
					" Exiting program." << endl;
				break;
			case 3:	cout << "At least one non-zero b_j is required."
					" Exiting program." << endl;
				break;
			case 4: cout << "'r' is the only character input that makes sense."
					" Exiting program." << endl;
				break;
			
			// Thrown in rdot() due to normalizing a zero vector
			case 5: cout << "The n vector must have non-zero magnitude."
				" Exiting program." << endl;
				break;
			
			// All other errors
			default:
				cout << "Indeterminate error. Algorithm not started." << endl;
		}
	}
	return 0;
}
