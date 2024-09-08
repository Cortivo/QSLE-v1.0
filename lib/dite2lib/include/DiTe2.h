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

#ifndef DITE2_H
#define DITE2_H

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
#include "inputParser.h"

class DiTe2 {

public:

	DiTe2(inputParser *ip0);
	virtual ~DiTe2();

	int calculateDiffusionTensor(std::string, std::string);
	double FrictionDihedral;

	mat getFrictionTensor(void);
	mat getDiffusionTensor(void);

private:

	mat c_alpha_per (vectorOfDoubles r);
	mat B (int number_rel_inter_coord, int* array, Col<int> vec_d, int **table, molecule mm);
	mat RP_sphere_diffusion (double Re, molecule mm);
	mat B_linear_generalization(int number_generalized, Col<int> number_internal, int **array, double **con_vec, Col<int> dih1, int **n_col, int ***table, molecule mm);

	inputParser *ip;	// Input parser object

	mat friction_tensor;	// Friction tensor matrix
	mat diffusion_tensor;	// Diffusion tensor matrix

};
#endif
