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

#include "DiTe2.h"


//Constants used in the code.

#define kb 1.38064852e-23	//Boltzmann constant expressed in J/K


/***************/
/* Constructor */
/***************/

DiTe2::DiTe2(inputParser *ip0)
{
	ip = ip0;
	return;
}


/*****************/
/* Deconstructor */
/*****************/

DiTe2::~DiTe2()
{
	return;
}


// Function that returns the matrix associated to the vector product:
// c_alpha_per w = w vec c_alpha

mat DiTe2::c_alpha_per (vectorOfDoubles r)
{
	mat ca(3,3);
	ca = { {0., r[2], -r[1]},
		   {-r[2], 0., r[0]},
		   {r[1], -r[0], 0.} };
	return ca;
}


// Matrix for the passage unconstrained-constrained (diffusion tensors)

mat DiTe2::B (int number_rel_inter_coord, int* array, Col<int> vec_d, int **table, molecule mm)
{
	int number_of_atoms = mm.getNAtoms();
	mat BB (3 * number_of_atoms, 6 + number_rel_inter_coord, fill::zeros);
	mat I_3(3, 3, fill::eye); // Identity matrix

	//WARNING: The atoms are numerated from 1
	for (int i = 0; i < 3 * number_of_atoms; i += 3) //N.B. Armadillo reads the matrix column by column. So I have to break the cycles.
	{
		//Traslational contribution

		BB.submat(i, 0, i + 2, 2) = I_3;


		//Rotational contribution

		BB.submat(i, 3, i + 2, 5) = c_alpha_per(mm.getAtomPosition(i / 3 + 1));
	}	

	int delay = 0;				

	for (int j = 0; j < number_rel_inter_coord; j++) //cycle on the columns
	{	
		if ((array[2 * j + 1] == 1) || (array[2 * j + 1] == 2)) //It holds only for bond lengths and bond angles
		{
			delay+=1; //number of both bond lenths and bond angles before the dihedral angles  	

			for (int i = 0; i < 3 * number_of_atoms; i += 3)
			{
				//Internal contribution

				vec a = mm.getFirstDerivative((i / 3 + 1), array[2 * j + 1], array[2 * j]);	
				BB(span(i, i + 2), j + 6) = a;
				a.reset();	
			}
		}
		else if ((array[2 * j + 1] == 3))
		{
			for (int i = 0; i < 3 * number_of_atoms; i += 3)
			{
				vec b(3, fill::zeros);
				for (int k = 0; k < vec_d(j - delay); k++) //cycle on the colums of the table
					{
						if ((table[j - delay][k] != 1) && (table[j - delay][k] != 2) && (table[j - delay][k] != 3))
						{
							b += mm.getFirstDerivative((i / 3 +1), 3, table[j - delay][k]);
						}
					}
					BB(span(i, i + 2), j + 6) = b;
					b.reset();
			}
		}
	}

	I_3.reset();		

	return BB;	
}

	
// Adimensional version of the diffusion tensor for the unconstrained spheres

mat DiTe2::RP_sphere_diffusion (double Re, molecule mm)
{
	int number_of_atoms = mm.getNAtoms();
	double mult;

	mat RP(3 * number_of_atoms, 3 * number_of_atoms, fill::zeros);
	mat I_3(3, 3, fill::eye); //Identity matrix

	for (int j = 0; j < 3 * number_of_atoms; j += 3) //cycle on the columns
	{
		RP.submat(j, j, j + 2, j + 2) = I_3;

		for (int i = 0; i < 3 * number_of_atoms; i += 3) //cycle on the raws
		{

			if (j !=i )
			{
				vec r_i = mm.getAtomPosition(i / 3 + 1);
				vec r_j = mm.getAtomPosition(j / 3 + 1);
				vec distance_ij = r_i - r_j;
				double n = norm(distance_ij); //Distance between i atom and j atom
				mat dyadic_ij = distance_ij * distance_ij.t();

				if (n > 2.0 * Re)
				{
					mult = (1. - 2. * ((Re * Re) / (n * n)));
					RP.submat(i, j, i + 2, j + 2) = (3. * Re / (4. * n * n * n)) * ((n * n + (2. / 3.) * Re * Re) * I_3 + mult * dyadic_ij);
				}
				else if (n <= 2. * Re)
				{
					RP.submat(i, j, i + 2, j + 2) = (1. - (9. / 32.) * (n/Re)) * I_3 + (3. / (32. * n * Re)) * dyadic_ij;	
				}
			}						
		}		
	}

	I_3.reset();

	return RP;
}

mat DiTe2::B_linear_generalization(int number_generalized, Col<int> number_internal, int **array, double **con_vec, Col<int> dih1, int **n_col, int ***table, molecule mm)
{

	int number_of_atoms = mm.getNAtoms();

	mat BB(3 * number_of_atoms, 6 + number_generalized, fill::zeros);
	mat I_3(3, 3, fill::eye); //Identity matrix

	
	//WARNING: The atoms are numerated from 1

	for (int i = 0; i < 3 * number_of_atoms; i += 3) // N.B. Armadillo reads the matrix column by column. So I have to break the cycles.
	{
		//Traslational contribution

		BB.submat(i, 0, i + 2, 2) = I_3;


		//Rotational contribution

		BB.submat(i, 3, i + 2, 5) = c_alpha_per(mm.getAtomPosition(i / 3 + 1));
	}	
		
	for (int i = 0; i < 3 * number_of_atoms; i += 3)
	{
		vec b(3, fill::zeros);
		int delay_r = 0;

		for (int l = 0; l < number_generalized; l++)
		{
			if (dih1(l) == 0)
			{
				delay_r += 1;
			}

			int delay = 0;

			for (int j = 0; j < number_internal(l) / 2; j++)
			{
				if ((array[l][2 * j + 1] == 1) || (array[l][2 * j + 1] == 2))
				{
					delay += 1;


					//Internal contribution for bond lengths and bond angles

					vec a = mm.getFirstDerivative((i / 3 + 1), array[l][2 * j + 1], array[l][2 * j]);
					a *= (1. / con_vec[l][j]);
					b += a;
					BB(span(i, i + 2), 6 + l) = b;
				}
				else if ((array[l][2 * j + 1] == 3))
				{
					vec a(3, fill::zeros);
					for (int k = 0; k < n_col[l - delay_r][j - delay]; k++) //cycle on the colums of the table
					{
						if ((table[l - delay_r][j - delay][k] !=1 ) && (table[l - delay_r][j - delay][k] !=2 ) && (table[l - delay_r][j - delay][k] != 3))
						{
							a += mm.getFirstDerivative((i / 3 + 1), 3, table[l - delay_r][j-delay][k]); // questo ritardo devo mettelo anche su un altro indice
						}
						a *= (1. / con_vec[l][j]);
					}

					b += a;
					BB(span(i, i + 2), 6 + l) = b;
				}
			}
		}
	}			
	
	I_3.reset();

	return BB;	
}	
	

/****************/
/* MAIN ROUTINE */	
/****************/

int DiTe2::calculateDiffusionTensor(std::string out1Name, std::string out2Name)
{
	// Run friction/diffusion tenrsors calculation

	if (ip->getQType() == 1)
	{
		//Construction of the B matrix
		#ifdef DEBUG
		std::cout << "Construction of the B matrix for the expression of the friction tensor\n";
		#endif
		//Friction and diffusion tensors

		mat d = RP_sphere_diffusion(ip->getReff(), ip->getMolecule());
#ifdef OUTPUT_SMALL_D
		std::cout << d << "\n\n";
#endif

		mat BB = B(ip->getNQ(), ip->getRelevantArray(), ip->getNumberCol(), ip->getDihedralTable() , ip->getMolecule());
#ifdef OUTPUT_BB
		std::cout << BB <<"\n \n";
#endif

		mat friction= BB.t() * d.i() * BB * (ip->getC() * ip->getViscosity() * 1.0e-10 * M_PI * ip->getReff());
		d.reset();
		friction.save(out1Name, arma_ascii);
		FrictionDihedral = friction(6,6);
		
		BB.reset();
		mat diffusion_tensor=(kb * ip->getTemperature() / (ip->getC() * ip->getViscosity() * 1.0e-10 * M_PI * ip->getReff())) * friction.i();
		friction.reset();
		diffusion_tensor.save(out2Name, arma_ascii);
		diffusion_tensor.reset();
	}
	else if (ip->getQType() == 2)
	{
		std::cout << "Framework of generalized coordinates\n";

		
		B_linear_generalization(ip->getNQ(), ip->getNInternal(), ip->getGeneralizedTable(), ip->getGeneralizedConstants(), ip->getNumberDihedral1(), ip->getNumberColCollective(), ip->getDihedralTableCollective(), ip->getMolecule());  


		//Friction and diffusion tensors


		mat d = RP_sphere_diffusion(ip->getReff(), ip->getMolecule());
#ifdef OUTPUT_SMALL_D
		std::cout << d << "\n\n";
#endif

		mat BB = B_linear_generalization(ip->getNQ(), ip->getNInternal(), ip->getGeneralizedTable(), ip->getGeneralizedConstants(), ip->getNumberDihedral1(), ip->getNumberColCollective(), ip->getDihedralTableCollective(), ip->getMolecule());
#ifdef OUTPUT_BB
		std::cout << BB <<"\n \n";
#endif

		friction_tensor = BB.t() * d.i() * BB * (ip->getC() * ip->getViscosity() * 1.0e-10 * M_PI * ip->getReff());
		d.reset();
		friction_tensor.save(out1Name, arma_ascii);
		friction_tensor = friction_tensor / (ip->getC() * ip->getViscosity() * 1.0e-10 * M_PI * ip->getReff());
	
		BB.reset();
		diffusion_tensor = (kb * ip->getTemperature() / (ip->getC() * ip->getViscosity() * 1.0e-10 * M_PI * ip->getReff())) * friction_tensor.i();
		friction_tensor.reset();
		diffusion_tensor.save(out2Name, arma_ascii);
		diffusion_tensor.reset();
	
	}
	else
	{
		std::cerr << "ERROR: Type of internal variable not recognized.\n";
		return 1;
	}
	std::cout<< std::endl;
	std::cout<< "### WARNING: DiTe2 standard output is suppressed to avoid verbosity" << std::endl << std::endl;
	/*
	std::cout << "=======================================================================" << std::endl;
	std::cout<<"DiTe2 results are stored in:" << std::endl;
	std::cout << "- " << out1Name << std::endl;
	std::cout << "- " << out2Name << std::endl ;
	std::cout << "=======================================================================" << std::endl<< std::endl;
	
	std::cout << "=======================================================================" << std::endl;
	std::cout << "For more information about DiTe2 program:" << std::endl;
	std::cout << "- https://doi.org/10.1002/jcc.25742"  << std::endl;
	std::cout << "- https://wwwdisc.chimica.unipd.it/mirco.zerbetto/?p=172"  << std::endl;
	std::cout << "=======================================================================" << std::endl << std::endl;
	*/
	
	return 0;
}


// Retrive data //

mat DiTe2::getFrictionTensor(void)
{
	return friction_tensor;
}

mat DiTe2::getDiffusionTensor(void)
{
	return diffusion_tensor;
}

