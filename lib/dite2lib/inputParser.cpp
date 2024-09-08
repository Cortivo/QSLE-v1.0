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

#include "inputParser.h"

/***************/
/* Constructor */
/***************/
inputParser::inputParser()
{
	// Some defaults: water at 298.15 K

	C = 6.0;
	Reff = 2.0;
	temperature = 298.15;
	viscosity = 8.94e-4;


	// Allocate arryay of IDs of 4 reference atoms

	ids = new int[4];


	// Init other quantities

	nUserBonds = 0;
	userBonds = NULL;

	pdbFileNameFound = 0;
	idsFound = 0;
	ReffFound = 0;
	CFound = 0;
	viscosityFound = 0;
	temperatureFound = 0;

	return;
}


/*****************/
/* Deconstructor */
/*****************/
inputParser::~inputParser()
{
	return;
}


// Parse input file
void inputParser::parseInputFile(std::string inputFileName)
{
	int tmpInt;
	int coordinatesFound = 0;

	std::fstream f;
	f.open(inputFileName, std::ios::in);
	std::string line, keyword;
	std::stringstream sline(std::ios_base::in);

	std::cout << "****************************" << std::endl;
	std::cout << "* Parsing DiTe2 input file *" << std::endl;
	std::cout << "****************************" << std::endl;

	while (getline(f, line))
	{
		sline.clear();
		sline.str(line);
		sline >> keyword;
		transform(keyword.begin(), keyword.end(), keyword.begin(), ::tolower);
		if(!keyword.compare(0, 1, "#"));	// skip the comment
		else if (!keyword.compare("pdb"))
		{
			sline >> pdbFileName;
			pdbFileNameFound = 1;
			std::cout << "- Input PDB file: " << pdbFileName << std::endl;
		}
		else if (!keyword.compare("refatoms"))
		{
			sline >> ids[0] >> ids[1] >> ids[2] >> ids[3];
			idsFound = 1;
			std::cout <<"- ID of reference atoms: " << ids[0] << " " << ids[1] << " " << ids[2] << " " << ids[3] << std::endl;
		}
		else if (!keyword.compare("reff"))
		{
			sline >> Reff;
			ReffFound = 1;
			std::cout << "- Effective radius for friction tensor: " << Reff << " Angstroms" << std::endl;
		}
		else if (!keyword.compare("c"))
		{
			sline >> C;
			CFound = 1;
			std::cout <<"- C (hydrodynamics boundary conditions): " << C << std::endl;
		}
		else if (!keyword.compare("viscosity"))
		{
			sline >> viscosity;
			viscosityFound = 1;
			std::cout << "- Medium viscosity: " << viscosity << " Pa s" << std::endl;
		}
		else if (!keyword.compare("temperature"))
		{
			sline >> temperature;
			temperatureFound = 1;
			std::cout << "- Temperature: " << temperature << " K" << std::endl;
		}
		else if (!keyword.compare("bonds"))
		{
			sline >> nUserBonds;
			std::cout << "Adding " << nUserBonds << " user defined bonds between atoms:" << std::endl;
			userBonds = new int[2 * nUserBonds];
			for (int i = 0; i < nUserBonds; ++i)
			{
				sline >> userBonds[i * 2 + 0] >> userBonds[i * 2 + 1];
				std::cout << "  " << userBonds[i * 2 + 0] << " and " << userBonds[i * 2 + 1] << std::endl;
			}
		}
		

	}
	qType = "simple";
	nq = 1;	
	coordinatesFound = 1;
	std::cout << "- Coordinates type: " << qType << std::endl;
	std::cout << "- Number of coordinates: " << nq << std::endl;

	relevant_array = new int[2 * nq];
	scan = new int[nq];
	firstLast = new double [2 * nq];
	relevant_array[0] = 4;
	relevant_array[1] = 3;
	scan[0] = 0;


	std::cout << "**********************" << std::endl;
	std::cout << "* End of DiTe2 input *" << std::endl;
	std::cout << "**********************" << std::endl;


	// Calculate the shifts for scan

	delta = new double[nq];
	for (int i = 0; i < nq; ++i)
	{
		if (scan[i] > 0)
			delta[i] = (firstLast[i * 2 + 1] - firstLast[i * 2 + 0]) / ((double)(scan[i]));
	}
	
	// Load the molecule

	loadMolecule(nUserBonds, userBonds);
	

	// Check consistency of all input information

	checkInputConsistency();
	if (!qType.compare("simple"))
	{
		coherence_internal();
		check_dihedrals_simple();
	}
	else if (!qType.compare("collective"))
	{
		coherence_generalized();
		check_dihedrals_collective();
	}	

	return;
}


// Check fundamental information

void inputParser::checkInputConsistency()
{
	int nerr = 0;
	if (!pdbFileNameFound)
	{
		std::cout << "ERROR: a PDB file must be specified via the ''pdb'' keyword in the input file." << std::endl;
		nerr ++;
	}
	if (!idsFound)
	{
		std::cout << "ERROR: the 4 reference atoms must be specified using the ''refAtoms ID1 ID2 ID3 ID4'' keyword." << std::endl;
		nerr ++;
	}
	if (!ReffFound)
	{
		std::cout << "WARNING: no Reff specified. Using default value of 2.0 Angstroms." << std::endl;
	}
	if (!CFound)
	{
		std::cout << "WARNING: no C specified. Uding default value of 6 (stick boundary conditions)." << std::endl;
	}
	if (!viscosityFound)
	{
		std::cout << "ERROR: viscosity (in Pa s) must be specified via the keyword ''viscosity''." << std::endl;
		nerr ++;
	}
	if (!temperatureFound)
	{
		std::cout << "WARNING: no temperature specified. Using default value of 298.15 K." << std::endl;
	}

	if (nerr)
	{
		std::cout << std::endl << "Terminating DiTe2 because of errors." << std::endl << std::endl;
		exit(1);
	}
	return;
}


// Load the molecule

void inputParser::loadMolecule(int nUserBonds, int *userBonds)
{
	// Configure where some useful information can be found

	config conf;
	std::string str;
        char* dite2Home;
        dite2Home = getenv ("QSLEINFO");
        if (dite2Home == NULL)
        {
                std::cout << std::endl << "### ERROR: the QSLEINFO envinronment variable is not set ###" << std::endl << std::endl;
                exit(1);
        }
        str.assign(dite2Home); str.append("/vdw.dat");
        conf.setVdwFile(str);
        str.assign(dite2Home); str.append("/atoms_masses.dat");
	conf.setMassFile(str);


	// Load the molecule

	mol = molecule(&conf);
	mol.setIOFileName(pdbFileName);
	mol.loadMolecule();


	// Grab some information

	int natoms = mol.getNAtoms();
	int nx = 3 * natoms;	//total number of coordinates
	int nIntTot = nx - 6;	//number of internal coordinates

	mol.setMainDihedralAngle(ids[0], ids[1], ids[2], ids[3]);


	// Build the Z-Matrix and prints it on the screen

	mol.buildZMatrix();


	// Express XYZ in AF

	mol.buildXYZfromZMatrix();


	// Dump the XYZ coordinates in AF

	mol.setIOFileName(pdbFileName + ".rebuild.xyz");
        mol.dumpMolecule();

	return;
}


// Check the consistency of the simple coordinates
	 
void inputParser::coherence_internal(void)
{
	int number_of_atoms = mol.getNAtoms();
	
	if (nq > 3 * number_of_atoms - 5 || nq < 0)
	{
		std::cerr << "ERROR: the number of coordinates must be N > 0 and N < 3*nAtoms - 5" << std::endl;
		exit(1);
	}
	

	// Consistency with the Z-matrix

	for (int i = 0; i < 2 * nq; i += 2)
	{
		if ((relevant_array[i] == 1) || ((relevant_array[i] == 2) && ((relevant_array[i + 1] == 2) || (relevant_array[i] == 3))) || ((relevant_array[i] == 3) && (relevant_array[i + 1] == 3)))
		{
			std::cout << "In the input file there are internal degrees not consistent with the Z-matrix" << std::endl;
			exit(1);
		}
	}


	// Check repetitions

	for (int i = 0; i < 2 * nq; i += 2)
	{
		for (int j = 2; j < 2 * nq; j += 2)
		{
			if ((i != j) && (relevant_array[i] == relevant_array[j]) && (relevant_array[i + 1] == relevant_array[j + 1]))
			{
				std::cout << "ERROR: the same internal degree is reapeated."<< std::endl;
				exit(1);
			}
		}
	}

	atom atom_i;
	zmatline z_i;
	atom atom_j;
	zmatline z_j;
	for (int i = 0; i < nq * 2; i += 2)
	{
		if (relevant_array[i + 1] == 3)
		{
			atom_i = mol.getAtom(relevant_array[i]);
			z_i = atom_i.getZMatrixEntry();
			{
				for (int j = 0; j < nq * 2; j += 2)
				{
					if (j != i)
					{	//This control works: proved with other couples of atoms
						atom_j = mol.getAtom(relevant_array[j]);
						z_j = atom_j.getZMatrixEntry();
						if ((z_j.A == z_i.A) && (z_j.B == z_i.B) && (relevant_array[j + 1] == 3))
						{
							std::cout << "The dihedral angles referenced by atoms\t" << relevant_array[i] << "\t and \t" << relevant_array[j] << "\t share the central atoms.\n";
							std::cout<<"It is not allowed to consider them as separate internal coordinate.\n";
							exit(1);
						}
					}
				}
			}
		}
	}

	return;
}


void inputParser::check_dihedrals_simple(void)
{
	number_dihedral = 0;
	vec dihedral1 = zeros<vec>(nq);

	for (int i = 0; i < nq * 2; i += 2)
	{
		if (relevant_array[i + 1] == 3)
		{
			number_dihedral ++; //number of the dihedral angles presented
			dihedral1(i >> 1) = relevant_array[i];
		}
	}

	int n_zeros = 0;
	Col<int> dihedral(number_dihedral, fill::zeros); // array of reference atom for the dihedral angles

	for (int i = 0; i < nq; ++i)
	{
		if (!dihedral1(i))
			n_zeros ++; //It counts the number of zeros
		else
			dihedral(i - n_zeros) = dihedral1(i);
	}

	dihedral1.reset();
	
	atom atom_i;
	zmatline z_i;
	atom atom_j;
	zmatline z_j;
	
	dihedral_table = new int*[number_dihedral];
	number_col = Col<int>(number_dihedral, fill::zeros);

	int col;

	for (int i = 0; i < number_dihedral; ++i)
	{
		atom_i = mol.getAtom(dihedral(i));
		z_i = atom_i.getZMatrixEntry();
		col = 1;
		dihedral_table[i] = new int[1];
		dihedral_table[i][0] = dihedral(i);
		number_col(i) = 1;

		for (int j = 1; j < mol.getNAtoms() + 1; ++j)
		{
			if (j != dihedral(i))
			{
				atom_j = mol.getAtom(j);
				z_j = atom_j.getZMatrixEntry();

				if ((z_j.A == z_i.A) && (z_j.B == z_i.B))
				{
					col ++;
					dihedral_table[i] = new int[col];
					dihedral_table[i][0] = dihedral(i);
					number_col(i) = col;

					for (int k = 1; k < col; k++)
						dihedral_table[i][k] = j;
				}
			}
		}
	}


	//Printing the control array for the dihedral angles 	
	//USE IT TO DEBUG	

	/*for(int i=0; i<number_dihedral; i++)    
	{
		for(int j=0; j<number_col(i); j++)
		{
			cout << dihedral_table[i][j]  << "  ";
		}
	cout << endl;
	}*/

	return;
}

// Check consistency for collective variables

void inputParser::coherence_generalized (void)
{
	int number_of_atoms = mol.getNAtoms();

	atom atom_i;
	zmatline z_i;
	atom atom_ii;
	zmatline z_ii;
	if ((nq > 3 * number_of_atoms - 5) || (nq < 0))
	{
		std::cerr << "ERROR: the number of coordinates must be N > 0 and N < 3*nAtoms - 5" << std::endl;
		exit(1);
	}
		
	for (int i = 0; i < nq; i++)
	{		
			if ((n_internal(i) > 3 * number_of_atoms - 5) || (n_internal(i) < 0))
			{
				std::cerr << "ERROR: the number of simple coordinates defining the collective coordinate must be between 1 and 3*nAtoms - 6" << std::endl;
				exit(1);
			}
	}

	for (int i = 0; i < nq; i++)
	{
		for (int j = 0; j < n_internal(i); j += 2)
		{
			if ((generalized_table[i][j] == 1) || ((generalized_table[i][j] == 2) && (generalized_table[i][j + 1] == 2)) || ((generalized_table[i][j] == 2) && (generalized_table[i][j + 1] == 3)) || ((generalized_table[i][j] == 3) && (generalized_table[i][j + 1] == 3)))
				{
					std::cout << "The choice of the internal degrees is not consistent. Please modify the input file\n";
					exit(1);
				}
		}
	}

	for (int i = 0; i < nq; i++)
	{
		for (int k = 0; k < n_internal(i); k += 2)
		{
			for (int m = 0; m < n_internal(i); m += 2)
			{
				if ((m != k) && (generalized_table[i][k] == generalized_table[i][m]) && (generalized_table[i][m + 1] == generalized_table[i][k + 1]))
				{
					std::cout<<"ERROR: a collective variable is defined twice.\n";
					exit(1);
				}
			}
		}
	}

	for (int i = 0; i < nq; i++)
	{
		for (int j = 0; j < nq; j++)
		{
			if (i != j)
			{
				for (int m = 0; m < n_internal(i); m += 2)
				{
					if ((generalized_table[i][m] == generalized_table[j][m]) && (generalized_table[i][m + 1] == generalized_table[j][m + 1]))
					{
					std::cout<<"ERROR: one simple coordinate is repeated twice in the definition of a collective coordinate.\n";
						exit(1);
					}
				}
			}
		}
	}

	for (int i = 0; i < nq; i++)
	{
		for (int j = 0; j < n_internal(i); j += 2)
		{
			if (generalized_table[i][j + 1] == 3)
			{
				atom_i=mol.getAtom(generalized_table[i][j]);
				z_i=atom_i.getZMatrixEntry();

				for (int ii = 0; ii < nq; ii++)
				{
					for (int jj = 0; jj < n_internal(ii); jj += 2)
					{
						if (((ii == i) && (jj != j)) || (ii != i) && (generalized_table[i][j] != generalized_table[ii][jj]))
						{
							atom_ii = mol.getAtom(generalized_table[ii][jj]);
							z_ii=atom_ii.getZMatrixEntry();

							if ((z_ii.A == z_i.A) && (z_ii.B == z_i.B) && (generalized_table[ii][jj + 1] == 3))
							{
								std::cout<<"The dihedral angles referenced by atoms\t"<<generalized_table[i][j]<<"\t and\t"<<generalized_table[ii][jj]<<"\t share the central atoms\n";
								std::cout<<"It is not allowed to consider them as separate internal coordinate.\n";
								exit(1);
							}
						}
					}
				}
			}
		}
	}
	
	return;
}

void inputParser::check_dihedrals_collective(void)
{
	std::cout << "Check the coherence of the dihedral angles chosen\n";
	number_dihedral1 = Col<int>(nq, fill::zeros);
		
	for (int i = 0; i < nq; ++i)
	{
		number_dihedral1(i) = 0;

		for (int m = 0; m < n_internal(i); m += 2)
		{
			if (generalized_table[i][m + 1] == 3)
				number_dihedral1(i) ++;
		}
	}
		
	int n_gen_dihedral = 0;	//number of generalized coordinates containing dihedral angles

	for (int i = 0; i < nq; i++)
	{
		if (number_dihedral1(i))
			n_gen_dihedral ++;
	}

	int n_zeros = 0;
	Col<int> number_dihedral(n_gen_dihedral, fill::zeros);	//for each generalized coordinate, how many dihedral angles are presented? Zeros are purged

	for (int i = 0; i < nq; i++)
	{
		if (!number_dihedral1(i))
			n_zeros ++;		//It counts the number of zeros
		else
			number_dihedral(i - n_zeros) = number_dihedral1(i);
	}

	//number_dihedral1.reset();
	
	int **dihedral = new int *[n_gen_dihedral];

	for (int i = 0; i < n_gen_dihedral; i++)
		dihedral[i] = new int[number_dihedral(i)];
	
	int l = 0, mm = 0;
	for (int i = 0; i < nq; i++)
	{
		if (!number_dihedral1(i))
			l ++;

		mm = 0;

		for (int j = 0; j < n_internal(i); j += 2)
		{
			if ((generalized_table[i][j + 1] == 1) || (generalized_table[i][j + 1] == 2))
				mm ++;
			else if (generalized_table[i][j + 1] == 3)
				dihedral[i - l][(j >> 1) - mm] = generalized_table[i][j];
		}
	}

	for(int i = 0; i < n_gen_dihedral; i++)
	{
		for(int j = 0; j < number_dihedral(i); j++)
		{
			std::cout << dihedral[i][j]  << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
		
	atom atom_ij;
	zmatline z_ij;
	atom atom_m;
	zmatline z_m;
	
	dihedral_table_coll = new int**[n_gen_dihedral];

	for (int i = 0; i < n_gen_dihedral; i++)
		dihedral_table_coll[i] = new int *[number_dihedral(i)];

	number_col_collective = new int*[n_gen_dihedral];

	for (int i = 0; i < n_gen_dihedral; ++i)
		number_col_collective[i] = new int[number_dihedral(i)];	
	
	int col;
	
	for (int i = 0; i < n_gen_dihedral; ++i)
	{
		for (int j = 0; j < number_dihedral(i); ++j)
		{
			atom_ij = mol.getAtom(dihedral[i][j]);
			z_ij = atom_ij.getZMatrixEntry();
			col = 1;
			dihedral_table_coll[i][j] = new int[1];
			dihedral_table_coll[i][j][0] = dihedral[i][j];
			number_col_collective[i][j] = 1;

			for (int m = 1 ; m < mol.getNAtoms() + 1; m++)	//cycle over atoms
			{
				if (m != dihedral[i][j])
				{
					atom_m = mol.getAtom(m);
					z_m = atom_m.getZMatrixEntry();
					if ((z_ij.A == z_m.A) && (z_ij.B == z_m.B))
					{
						col ++;
						dihedral_table_coll[i][j] = new int[col];
						dihedral_table_coll[i][j][0] = dihedral[i][j];
						number_col_collective[i][j] = col;
						for (int k = 1; k < col; k++)
							dihedral_table_coll[i][j][k] = m;
					}
				}
			}
		}
	}

	for(int i = 0; i < n_gen_dihedral; i++)    
	{
		for(int j = 0; j < number_dihedral(i); j++)
		{
			for (int k = 0; k < number_col_collective[i][j]; k++)
				std::cout << dihedral_table_coll[i][j][k] << "  ";
			std::cout << std::endl;
		}
			std::cout << std::endl;
	}

	return;
}


// Methods to retrive data

molecule inputParser::getMolecule(void)
{
	return mol;
}

int inputParser::getQType(void)
{
	if (!qType.compare("simple"))
		return 1;
	else if (!qType.compare("collective"))
		return 2;
	else
		return 0;
}

int inputParser::getNQ(void)
{
	return nq;
}

int* inputParser::getRelevantArray(void)
{
	return relevant_array;
}

int inputParser::getRelevantArray(int i)
{
	return relevant_array[i];
}

int** inputParser::getDihedralTable(void)
{
	return dihedral_table;
}

Col<int> inputParser::getNumberCol(void)
{
	return number_col;
}

int** inputParser::getGeneralizedTable(void)
{
	return generalized_table;
}

int inputParser::getGeneralizedTable(int i, int j)
{
	return generalized_table[i][j];
}

double** inputParser::getGeneralizedConstants(void)
{
	return generalized_constants;
}

double inputParser::getGeneralizedConstants(int i, int j)
{
	return generalized_constants[i][j];
}

Col<int> inputParser::getNumberDihedral1(void)
{
	return number_dihedral1;
}

int** inputParser::getNumberColCollective(void)
{
	return number_col_collective;
}

int*** inputParser::getDihedralTableCollective(void)
{
	return dihedral_table_coll;
}

Col<int> inputParser::getNInternal(void)
{
	return n_internal;
}

int inputParser::getNInternal(int i)
{
	return n_internal(i);
}

double inputParser::getReff(void)
{
	return Reff;
}

double inputParser::getC(void)
{
	return C;
}

double inputParser::getTemperature(void)
{
	return temperature;
}

double inputParser::getViscosity(void)
{
	return viscosity;
}

std::string inputParser::getPDBFileName(void)
{
	return pdbFileName;
}

int* inputParser::getScan(void)
{
	return scan;
}

int inputParser::getScan(int i)
{
	return scan[i];
}

double inputParser::getFirst(int i)
{
	return firstLast[i * 2 + 0];
}

double inputParser::getDelta(int i)
{
	return delta[i];
}


// Methods to act on the molecule

void inputParser::changeZMatrix(int i, double val)
{
	if (!qType.compare("simple"))
	{
		int aID = relevant_array[i * 2 + 0];
		int ABC = relevant_array[i * 2 + 1];
		zmatline z = mol.getAtom(aID).getZMatrixEntry();
		switch (ABC)
		{
			case 1:
			{
				z.d = val;
				mol.changeZMatrix(aID, z);
				break;
			}
			case 2:
			{
				z.theta = val * M_PI / 180.0;
				mol.changeZMatrix(aID, z);
				break;
			}
			default: // A.K.A. case 3
			{
				val *= M_PI / 180.0;

				double deltaPhi = val - z.phi;

				z.phi = val;
				mol.changeZMatrix(aID, z);
	

				// Rotate all the dihedrals sharing the same pair of central atoms

				int row = 0;
				for (int ir = 0; ir < number_dihedral; ++ir)
				{
					if (dihedral_table[ir][0] == aID)
					{
						row = ir;
						break;
					}
				}

				for (int j = 1; j < number_col[row]; ++j)
				{
					aID = dihedral_table[row][j];
					z = mol.getAtom(aID).getZMatrixEntry();
					z.phi += deltaPhi;
					mol.changeZMatrix(aID, z);
				}
				break;
			}
		}
		mol.outputZMatrix("screen");
	}
	else
	{
		int aID, ABC;
		zmatline z;
		for (int ni = 0; ni < n_internal(i); ni += 2)
		{
			aID = generalized_table[i][ni];
			ABC = generalized_table[i][ni + 1];
			z = mol.getAtom(aID).getZMatrixEntry();
			switch (ABC)
			{
				case 1:
				{
					z.d += delta[i] / generalized_constants[i][ni >> 1];
					mol.changeZMatrix(aID, z);
					break;
				}
				case 2:
				{
					z.theta += delta[i] / generalized_constants[i][ni >> 1];
					mol.changeZMatrix(aID, z);
					break;
				}
				default: // A.K.A. case 3
				{
					z.phi += delta[i] / generalized_constants[i][ni >> 1];
					mol.changeZMatrix(aID, z);


					// Rotate the torsion angles related to this one

					int row = 0;
					for (int ir = 0; ir < number_dihedral1(i); ++ir)
					{
						if (dihedral_table_coll[i][ir][0] == aID)
						{
							row = ir;
							break;
						}
					}
	
					for (int j = 1; j < number_col_collective[i][row]; ++j)
					{
						aID = dihedral_table_coll[i][row][j];
						z = mol.getAtom(aID).getZMatrixEntry();
						z.phi += delta[i] / generalized_constants[i][ni >> 1];
						mol.changeZMatrix(aID, z);
					}
					
					break;
				}
			}
		}
	}

	return;
}

void inputParser::dumpMolecule(std::string name)
{
	mol.setIOFileName(name);
	mol.dumpMolecule();
	return;
}

