#ifndef ATOM_H_
#define ATOM_H_

#include <cstdlib>
#include <ctype.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>

#include "types.h"

#include "matrix4D.h"

class atom
{
	public:	

		int id;  // Actual atom ID
		int id0; // Atom ID in the original XYZ or PDB

		int activeID;
		int activeID0;

		int getID(void);
		int getID0(void);

		atom();
		virtual ~atom();
	
		double x();
		void x(double);
		double y();
		void y(double);
		double z();
		void z(double);
		double x0();
		void x0(double);
		double y0();
		void y0(double);
		double z0();
		void z0(double);
		int isActive(void);
		void isActive(int);
		void savePosition(void);
		vectorOfDoubles getPosition();
		vectorOfDoubles getPosition0();
		void setPosition(vectorOfDoubles);
	
		int getNBonds();
		void setNBonds(int);
		int getBond(int);
		vectorOfIntegers getBonds();
		void setBonds(vectorOfIntegers);
		void addBond(int);
		void clearBonds();
	
		void setResidueNumber(int);
		int getResidueNumber();
		void setResidueName(std::string);
		std::string getResidueName();
		void setAtomType(std::string);
		std::string getAtomType();
		void setAtomTypePDB(std::string);
		std::string getAtomTypePDB();

		void setFrictionMultiplier(double);
		double getFrictionMutiplier(void);
	
		void setCellNumbers(int,int,int);
		void setCellNumbers(vectorOfIntegers);
		vectorOfIntegers getCellNumbers();

		void setMass(vectorOfStrings, vectorOfDoubles);
		double getMass();

		// Z-Matrix related stuff

		void setZMatrixEntry(int, double, int, double, int, double);
		void setZMatrixEntry(zmatline);
		matrix4D getBmat(void);
		void overrideBmat(zmatline);
		matrix4D getBmatFirstDerivative(int);
		matrix4D getBmatSecondDerivative(int, int);
		zmatline getZMatrixEntry();
		int getZMatAtom(int);
		double getZMatValue(int);
		void setAmat(matrix4D);
		matrix4D getAmat(void);
		void setChain(vectorOfIntegers);
		vectorOfIntegers getChain(void);
		int getChainSize(void);

		// Fragments

		void setNFrag(int);
		int getNFrag(void);
	
	private:
		int nBonds;
		int residueNumber;
		vectorOfIntegers cellNumbers;
		std::string atomType;
		std::string atomTypePDB;
		std::string residueName;
		double mass;
		vectorOfIntegers bonds;
		vectorOfDoubles position;
		vectorOfDoubles position0;

		// Z-Matrix related stuff
 
		vectorOfIntegers chain; // chain of atoms through reference atom (permuted order numbering), 0-based
		zmatline zmat;
		matrix4D Bmat; // Roto-translation matrix B
		matrix4D Bmat_d, Bmat_t, Bmat_p; // First derivatives of B
		matrix4D Bmat_dt, Bmat_dp, Bmat_td, Bmat_tt, Bmat_tp, Bmat_pd, Bmat_pt, Bmat_pp; // Second derivatives of B
		matrix4D Amat; // Tranformation matrix from actual atom to reference atom
		void buildB(void);
		void buildB(zmatline);
 
		// Fragmentation

		int nfrag, next;

		// Active flag:
		// 0 = the atom is excluded from the model
		// 1 = the atom is part of the model
		int active;

		// Friction multiplier
		double csiMult;
};

#endif
