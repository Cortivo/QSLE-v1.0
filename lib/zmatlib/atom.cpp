/*
 *  atom.cpp
 *  gc
 *
 *  Created by Mirco Zerbetto on 6/24/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "atom.h"

/***************/
/* Constructor */
/***************/
atom::atom()
{
	position.push_back(0.0);
	position.push_back(0.0);
	position.push_back(0.0);
	nBonds = 0;
	bonds = std::vector<int>(0);
	position0.push_back(0.0);
	position0.push_back(0.0);
	position0.push_back(0.0);

	residueNumber = 0;
	active        = 1;

	csiMult = 1.0;
}

/*****************/
/* Deconstructor */
/*****************/
atom::~atom()
{
}

/***********************/
/* Method to return ID */
/***********************/

int atom::getID(void)
{
	return id;
}

int atom::getID0(void)
{
	return id0;
}

/*******************************/
/* Return x, y or z coordinate */
/*******************************/
double atom::x()
{
	return position[0];
}
double atom::y()
{
	return position[1];
}
double atom::z()
{
	return position[2];
}

double atom::x0()
{
	return position0[0];
}
double atom::y0()
{
	return position0[1];
}
double atom::z0()
{
	return position0[2];
}

/*******************************/
/* Set x, y or z coordinate */
/*******************************/
void atom::x(double val)
{
	position[0] = val;
	return;
}
void atom::y(double val)
{
	position[1] = val;
	return;
}
void atom::z(double val)
{
	position[2] = val;
	return;
}

void atom::x0(double val)
{
	position0[0] = val;
	return;
}
void atom::y0(double val)
{
	position0[1] = val;
	return;
}
void atom::z0(double val)
{
	position0[2] = val;
	return;
}

void atom::savePosition(void)
{
	position0[0] = position[0];
	position0[1] = position[1];
	position0[2] = position[2];
	return;
}

/*******************************/
/* Set/Get the ''active'' flag */
/*******************************/

void atom::isActive(int d)
{
	active = d;
	return;
}

int atom::isActive(void)
{
	return active;
}

/*******************/
/* Return position */
/*******************/
vectorOfDoubles atom::getPosition()
{
	return position;
}

vectorOfDoubles atom::getPosition0()
{
	return position0;
}

/****************/
/* Set position */
/****************/
void atom::setPosition(vectorOfDoubles p)
{
	if (p.size() < 1)
		std::cout << "WARNING >> atom::setPosition(vectorOfDoubles) : passed vector has size < 1. Position remains unchanged" << std::endl;
	else if (p.size() == 1)
	{
		std::cout << "WARNING >> atom::setPosition(vectorOfDoubles) : passed vector has size 1. Only x coordinate will be changed" << std::endl;
		position[0] = p[0];
	}
	else if (p.size() == 2)
	{
		std::cout << "WARNING >> atom::setPosition(vectorOfDoubles) : passed vector has size 2. Only x and y coordinates will be changed" << std::endl;
		position[0] = p[0];
		position[1] = p[1];
	}
	else if (p.size() == 3)
	{
		position[0] = p[0];
		position[1] = p[1];
		position[2] = p[2];
	}
	else
	{
		std::cout << "WARNING >> atom::setPosition(vectorOfDoubles) : passed vector has size > 3. Only first three elements will be used" << std::endl;
		position[0] = p[0];
		position[1] = p[1];
		position[2] = p[2];
	}
	
	return;
}

/******************************/
/* Return the number of bonds */
/******************************/
int atom::getNBonds()
{
	return nBonds;
}

/***************************/
/* Set the number of bonds */
/***************************/
void atom::setNBonds(int n)
{
	nBonds = n;
	return;
}

/**********************/
/* Return bonds array */
/**********************/
vectorOfIntegers atom::getBonds()
{
	return bonds;
}

int atom::getBond(int b)
{
	// b is 0-based index
	// return is 1-based ID
	return bonds.at(b);
}

/*******************************************************/
/* Set the bonds array (update nBonds for consistency) */
/*******************************************************/
void atom::setBonds(vectorOfIntegers b)
{
	int n = b.size();
	if (n != nBonds)
		std::cout << "WARNING >> atom::setBonds(vectorOfIntegers) : actual nBonds (" << nBonds << ") has been changed to passed array size (" << n << ") for consistency" << std::endl;
	nBonds = n;
	bonds = b;
}

/**************/
/* Add a bond */
/**************/
void atom::addBond(int at)
{
	nBonds++;
	bonds.push_back(at);
	return;
}

/***************/
/* Clear bonds */
/***************/
void atom::clearBonds()
{
	nBonds = 0;
	bonds.clear();
}

/******************/
/* Residue number */
/******************/
void atom::setResidueNumber(int n)
{
	residueNumber = n;
	return;
}

int atom::getResidueNumber()
{
	return residueNumber;
}

/****************/
/* Residue name */
/****************/
void atom::setResidueName(std::string r)
{
	residueName = r;
	return;
}

std::string atom::getResidueName()
{
	return residueName;
}

/*************/
/* Atom type */
/*************/
void atom::setAtomType(std::string a)
{
	atomType = a;
	return;
}
void atom::setAtomTypePDB(std::string a)
{
	atomTypePDB = a;
	atomType = a;
	int pos = 0;
	while (isdigit(atomType[pos]))
		pos++;
	pos++;

	if (atomType[pos] >= 'a' && atomType[pos] <= 'z')
		atomType = atomType.substr(pos-1, 2);
	else
		atomType = atomType.substr(pos-1, 1);
}

std::string atom::getAtomType()
{
	return atomType;
}

std::string atom::getAtomTypePDB()
{
	return atomTypePDB;
}

/***************/
/* Cell number */
/***************/

void atom::setCellNumbers(int c1, int c2, int c3)
{
	cellNumbers = vectorOfIntegers(3,0.0);
	cellNumbers[0] = c1;
	cellNumbers[1] = c2;
	cellNumbers[2] = c3;
	return;
}

void atom::setCellNumbers(vectorOfIntegers c)
{
	cellNumbers = c;
	return;
}

vectorOfIntegers atom::getCellNumbers()
{
	return cellNumbers;
}

/***************/
/* Atomic mass */
/***************/
void atom::setMass(vectorOfStrings types, vectorOfDoubles masses)
{
	unsigned int i = 0;
	for (i = 0; i < types.size(); i++)
		if (!atomType.compare(types[i]))
		{
			mass = masses[i];
			return;
		}
	std::cout << "ERROR : atom::setMass(vectorOfStrings, vectorOfDoubles) >> Cannot find mass for atom " << atomType << std::endl;
	exit(1);
}

double atom::getMass()
{
	return mass;
}

/************/
/* Z-Matrix */
/************/

void atom::setZMatrixEntry(int a, double d, int b, double t, int c, double p)
{
	zmat.A = a;
	zmat.B = b;
	zmat.C = c;
	zmat.d = d;
	zmat.theta = t;
	zmat.phi = p;
	buildB();
	return;
}

void atom::setZMatrixEntry(zmatline z)
{
	zmat = z;
	buildB();
	return;
}

matrix4D atom::getBmat(void)
{
	return Bmat;
}

matrix4D atom::getBmatFirstDerivative(int i)
{
	switch (i)
	{
		case 1:
		{
			return Bmat_d;
		}
		case 2:
		{
			return Bmat_t;
		}
		case 3:
		{
			return Bmat_p;
		}
		default:
		{
			std::cout << "ERROR : atom::getBmatFirstDerivative(int) >> index " << i << " out of range. Possible values:\n1 -> derivative with respect to bond length\n2 -> derivative with respect to bond angle\n3 -> derivative with respect to dihedral angle.\n\n" << std::endl;
			exit(1);
		}
	}
	return Bmat;
}


matrix4D atom::getBmatSecondDerivative(int i, int j)
{
	switch (i)
	{
		case 1:
		{
			switch (j)
			{
				case 1:
				{
					return (matrix4D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
				}
				case 2:
				{
					return Bmat_dt;
				}
				case 3:
				{
					return Bmat_dp;
				}
				default:
				{
					std::cout << "ERROR : atom::getBmatSecondDerivative(int, int) >> index of second differentiation," << j << ", out of range. Possible values:\n1 -> derivative with respect to bond length\n2 -> derivative with respect to bond angle\n3 -> derivative with respect to dihedral angle.\n\n" << std::endl;
			exit(1);
				}
			}
		}
		case 2:
		{
			switch (j)
			{
				case 1:
				{
					return Bmat_td;
				}
				case 2:
				{
					return Bmat_tt;
				}
				case 3:
				{
					return Bmat_tp;
				}
				default:
				{
					std::cout << "ERROR : atom::getBmatSecondDerivative(int, int) >> index of second differentiation," << j << ", out of range. Possible values:\n1 -> derivative with respect to bond length\n2 -> derivative with respect to bond angle\n3 -> derivative with respect to dihedral angle.\n\n" << std::endl;
			exit(1);
				}
			}
		}
		case 3:
		{
			switch (j)
			{
				case 1:
				{
					return Bmat_pd;
				}
				case 2:
				{
					return Bmat_pt;
				}
				case 3:
				{
					return Bmat_pp;
				}
				default:
				{
					std::cout << "ERROR : atom::getBmatSecondDerivative(int, int) >> index of second differentiation," << j << ", out of range. Possible values:\n1 -> derivative with respect to bond length\n2 -> derivative with respect to bond angle\n3 -> derivative with respect to dihedral angle.\n\n" << std::endl;
			exit(1);
				}
			}
		}
		default:
		{
			std::cout << "ERROR : atom::getBmatSecondDerivative(int, int) >> index of first differentiation," << i << ", out of range. Possible values:\n1 -> derivative with respect to bond length\n2 -> derivative with respect to bond angle\n3 -> derivative with respect to dihedral angle.\n\n" << std::endl;
			exit(1);
		}
	}

	return Bmat_d;
}

zmatline atom::getZMatrixEntry()
{
	return zmat;
}

int atom::getZMatAtom(int i)
{
	switch (i)
	{
		case 1:
			return zmat.A;
		case 2:
			return zmat.B;
		case 3:
			return zmat.C;
		default:
		{
			std::cout << "ERROR : atom::getZMAtAtom(int) >> index " << i << " can be 1, 2, or 3" << std::endl;
			exit(1);

		}
	}
	return 0;
}

double atom::getZMatValue(int i)
{
	switch (i)
	{
		case 1:
			return zmat.d;
		case 2:
			return zmat.theta;
		case 3:
			return zmat.phi;
		default:
		{
			std::cout << "ERROR : atom::getZMAtAtom(int) >> index " << i << " can be 1, 2, or 3" << std::endl;
			exit(1);

		}
	}
	return 0;
}

/* Methods to built the roto-translation matrix B and its
   first and second derivatives */

void atom::overrideBmat(zmatline overZ)
{
	buildB(overZ);
	return;
}

void atom::buildB(zmatline overZ)
{
	double ct, st, cp, sp;

	ct = cos(overZ.theta);
	st = sin(overZ.theta);
	cp = cos(overZ.phi);
	sp = sin(overZ.phi);

	// B matrix

	Bmat = matrix4D(-ct    , -st     ,  0.0, -overZ.d * ct     ,
			 st * cp, -ct * cp, -sp , overZ.d * st * cp,
			 st * sp, -ct * sp,  cp , overZ.d * st * sp,
			 0.0    ,  0.0    ,  0.0, 1.0              );

	// First derivatives

	Bmat_d = matrix4D( 0.0,  0.0,  0.0, -ct    ,
			   0.0,  0.0,  0.0,  st * cp,
			   0.0,  0.0,  0.0,  st * sp,
			   0.0,  0.0,  0.0,  0.0     );


	Bmat_t = matrix4D( st     , -ct     ,  0.0,  overZ.d * st     ,
			   ct * cp,  st * cp,  0.0,  overZ.d * ct * cp,
			   ct * sp,  st * sp,  0.0,  overZ.d * ct * sp,
			   0.0    ,  0.0    ,  0.0,  0.0              );


	Bmat_p = matrix4D( 0.0    ,  0.0    ,  0.0,  0.0             ,
			  -st * sp,  ct * sp, -cp , -overZ.d * st * sp,
			   st * cp, -ct * cp, -sp ,  overZ.d * st * cp,
			   0.0    ,  0.0    ,  0.0,  0.0              );

	// Second derivatives: order of differentiation Bmat_FirstSecond

	Bmat_dt = matrix4D( 0.0,  0.0,  0.0,  st    ,
			    0.0,  0.0,  0.0,  ct * cp,
			    0.0,  0.0,  0.0,  ct * sp,
			    0.0,  0.0,  0.0,  0.0     );

	Bmat_dp = matrix4D( 0.0,  0.0,  0.0,  0.0    ,
			    0.0,  0.0,  0.0, -st * sp,
			    0.0,  0.0,  0.0,  st * cp,
			    0.0,  0.0,  0.0,  0.0     );


	Bmat_td = matrix4D( 0.0,  0.0,  0.0,  st     ,
			    0.0,  0.0,  0.0,  ct * cp,
			    0.0,  0.0,  0.0,  ct * sp,
			    0.0,  0.0,  0.0,  0.0     );

	Bmat_tt = matrix4D( ct     ,  st     ,  0.0,  overZ.d * ct     ,
			   -st * cp,  ct * cp,  0.0, -overZ.d * st * cp,
			   -st * sp,  ct * sp,  0.0, -overZ.d * st * sp,
			    0.0    ,  0.0    ,  0.0,  0.0              );

	Bmat_tp = matrix4D( 0.0    ,  0.0    ,  0.0,  0.0     ,
			   -ct * sp, -st * sp,  0.0, -overZ.d * ct * sp,
			    ct * cp,  st * cp,  0.0,  overZ.d * ct * cp,
			    0.0    ,  0.0    ,  0.0,  0.0              );


	Bmat_pd = matrix4D( 0.0,  0.0,  0.0,  0.0    ,
			    0.0,  0.0,  0.0, -st * sp,
			    0.0,  0.0,  0.0,  st * cp,
			    0.0,  0.0,  0.0,  0.0     );

	Bmat_pt = matrix4D( 0.0    ,  0.0    ,  0.0,  0.0             ,
			   -ct * sp, -st * sp,  0.0, -overZ.d * ct * sp,
			    ct * cp,  st * cp,  0.0,  overZ.d * ct * cp,
			    0.0    ,  0.0    ,  0.0,  0.0              );

	Bmat_pp = matrix4D( 0.0    ,  0.0    ,  0.0,  0.0             ,
			   -st * cp,  ct * cp,  sp , -overZ.d * st * cp,
			   -st * sp,  ct * sp, -cp , -overZ.d * st * sp,
			    0.0    ,  0.0    ,  0.0,  0.0              );

	return;
}

void atom::buildB(void)
{
	double ct, st, cp, sp;

	ct = cos(zmat.theta);
	st = sin(zmat.theta);
	cp = cos(zmat.phi);
	sp = sin(zmat.phi);

	// B matrix

	Bmat = matrix4D(-ct    , -st     ,  0.0, -zmat.d * ct     ,
			 st * cp, -ct * cp, -sp , zmat.d * st * cp,
			 st * sp, -ct * sp,  cp , zmat.d * st * sp,
			 0.0    ,  0.0    ,  0.0, 1.0              );

	// First derivatives

	Bmat_d = matrix4D( 0.0,  0.0,  0.0, -ct    ,
			   0.0,  0.0,  0.0,  st * cp,
			   0.0,  0.0,  0.0,  st * sp,
			   0.0,  0.0,  0.0,  0.0     );


	Bmat_t = matrix4D( st     , -ct     ,  0.0,  zmat.d * st     ,
			   ct * cp,  st * cp,  0.0,  zmat.d * ct * cp,
			   ct * sp,  st * sp,  0.0,  zmat.d * ct * sp,
			   0.0    ,  0.0    ,  0.0,  0.0              );


	Bmat_p = matrix4D( 0.0    ,  0.0    ,  0.0,  0.0             ,
			  -st * sp,  ct * sp, -cp , -zmat.d * st * sp,
			   st * cp, -ct * cp, -sp ,  zmat.d * st * cp,
			   0.0    ,  0.0    ,  0.0,  0.0              );

	// Second derivatives: order of differentiation Bmat_FirstSecond

	Bmat_dt = matrix4D( 0.0,  0.0,  0.0,  st    ,
			    0.0,  0.0,  0.0,  ct * cp,
			    0.0,  0.0,  0.0,  ct * sp,
			    0.0,  0.0,  0.0,  0.0     );

	Bmat_dp = matrix4D( 0.0,  0.0,  0.0,  0.0    ,
			    0.0,  0.0,  0.0, -st * sp,
			    0.0,  0.0,  0.0,  st * cp,
			    0.0,  0.0,  0.0,  0.0     );


	Bmat_td = matrix4D( 0.0,  0.0,  0.0,  st     ,
			    0.0,  0.0,  0.0,  ct * cp,
			    0.0,  0.0,  0.0,  ct * sp,
			    0.0,  0.0,  0.0,  0.0     );

	Bmat_tt = matrix4D( ct     ,  st     ,  0.0,  zmat.d * ct     ,
			   -st * cp,  ct * cp,  0.0, -zmat.d * st * cp,
			   -st * sp,  ct * sp,  0.0, -zmat.d * st * sp,
			    0.0    ,  0.0    ,  0.0,  0.0              );

	Bmat_tp = matrix4D( 0.0    ,  0.0    ,  0.0,  0.0     ,
			   -ct * sp, -st * sp,  0.0, -zmat.d * ct * sp,
			    ct * cp,  st * cp,  0.0,  zmat.d * ct * cp,
			    0.0    ,  0.0    ,  0.0,  0.0              );


	Bmat_pd = matrix4D( 0.0,  0.0,  0.0,  0.0    ,
			    0.0,  0.0,  0.0, -st * sp,
			    0.0,  0.0,  0.0,  st * cp,
			    0.0,  0.0,  0.0,  0.0     );

	Bmat_pt = matrix4D( 0.0    ,  0.0    ,  0.0,  0.0             ,
			   -ct * sp, -st * sp,  0.0, -zmat.d * ct * sp,
			    ct * cp,  st * cp,  0.0,  zmat.d * ct * cp,
			    0.0    ,  0.0    ,  0.0,  0.0              );

	Bmat_pp = matrix4D( 0.0    ,  0.0    ,  0.0,  0.0             ,
			   -st * cp,  ct * cp,  sp , -zmat.d * st * cp,
			   -st * sp,  ct * sp, -cp , -zmat.d * st * sp,
			    0.0    ,  0.0    ,  0.0,  0.0              );

	return;
}

void atom::setAmat(matrix4D a)
{
	Amat = a;
	return;
}

matrix4D atom::getAmat(void)
{
	return Amat;
}

void atom::setChain(vectorOfIntegers c)
{
	chain = c;
	return;
}

vectorOfIntegers atom::getChain(void)
{
	return chain;
}

int atom::getChainSize(void)
{
	return chain.size();
}

/*****************/
/* Fragmnet data */
/*****************/

void atom::setNFrag(int n)
{
	nfrag = n;
	return;
}

int atom::getNFrag(void)
{
	return nfrag;
}

/***********************/
/* Friction multiplier */
/***********************/

void atom::setFrictionMultiplier(double m)
{
	csiMult = m;
	return;
}

double atom::getFrictionMutiplier(void)
{
	return csiMult;
}

