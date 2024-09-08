/***********************************************************************************
 * Dite 2 v1.0 - Program to evaluae full diffusion tensor of flexible molecules    *
 * Copyright (C) 2018  Jonathan Campeggio, Antonino Polimeno, Mirco Zerbetto       * 
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or any later version.                                           *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 ***********************************************************************************
 * Authors: Jonathan Campeggio, Antonino Polimeno, Mirco Zerbetto                  *
 * Dipartimento di Scienze Chimiche - Universita' di Padova - Italy                *
 * E-mail: mirco.zerbetto@unipd.it                                                 *
 ***********************************************************************************/

#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <array> 	//Module for the manipulation of fixed array. 
#include <vector>	//Module for the manipulation of dynamic array.

#include "armadillo"
using namespace arma;

#include "zmat.h"

class inputParser {

public:

	inputParser();
	virtual ~inputParser();

	void parseInputFile(std::string inputFileName);

	molecule getMolecule(void);

	int getQType(void);
	int getNQ(void);
	int* getRelevantArray(void);
	int getRelevantArray(int);
	int** getDihedralTable(void);
	Col<int>getNumberCol(void);
	int** getGeneralizedTable(void);
	int getGeneralizedTable(int i, int j);
	double** getGeneralizedConstants(void);
	double getGeneralizedConstants(int i, int j);
	Col<int>getNumberDihedral1(void);
	int** getNumberColCollective(void);
	int*** getDihedralTableCollective(void);
	Col<int> getNInternal(void);
	int getNInternal(int i);
	int* getScan(void);
	int getScan(int);
	double getFirst(int i);
	double getDelta(int i);

	double getReff(void);
	double getC(void);
	double getTemperature(void);
	double getViscosity(void);
	std::string getPDBFileName(void);

	void changeZMatrix(int i, double val);
	void dumpMolecule(std::string name);

private:

	void checkInputConsistency(void);
	void loadMolecule(int nUserBonds, int *userBonds);
	void coherence_internal(void);
	void check_dihedrals_simple(void);
	void coherence_generalized(void);
	void check_dihedrals_collective(void);

	int *ids;			// Array of IDs of 4 reference atoms
	int nq;				// Number of coordinates
	int *relevant_array;		// Array 2xN with the information on N simple coordinates
	int *scan;			// Number of discretization points for a scan of the coordinates
	int **generalized_table;	// Array with the definition of the collective variables
	Col<int>n_internal;		// Array of internal coordinates per collective coordinate
	int nUserBonds;			// Number of user-defined bonds
	int *userBonds;			// Array with user-defined pairs of atoms to be bonded
	int number_dihedral;		// Number of dihedral angles, simple coordinates
	int **dihedral_table;		// Table of dihedrals for simple coordinates
	int **number_col_collective;	//
	int ***dihedral_table_coll;	// Table of dihedrals for simple coordinates
	Col<int>number_col;		// Array with atoms indeces of dihedral angles
	Col<int>number_dihedral1;	// Number of dihedrals per collective coordinate 

	int pdbFileNameFound;
	int idsFound;
	int ReffFound;
	int CFound;
	int viscosityFound;
	int temperatureFound;

	double C; 			// hydrodynamic boundary conditions
	double Reff;			// Effective beads radius
	double temperature;		// Temperature
	double viscosity;		// Viscosity
	double *firstLast;		// First and last points for coordinate scan
	double *delta;			// Array with variations of coordinates during scan
	double **generalized_constants;	// Array with the coefficients of the collective variables

	std::string qType;		// Coordinates type
	std::string pdbFileName;	// Input PDB file

	molecule mol;			// Molecule

};
#endif
