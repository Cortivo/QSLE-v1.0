#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

#include "atom.h"
#include "config.h"
#include "frameRT.h"
#include "jacobi.h"

#include "types.h"
#include "globals.h"

typedef std::vector<atom> vectorOfAtoms;

class molecule
{
	
public:
	molecule();
	molecule(config *);
	virtual ~molecule();
	void setIOFileName(std::string);
	std::string getFileName();
	int loadMolecule();
	int dumpMolecule();
	
	std::string toString(bool,bool,bool,bool,bool,bool);
	
	int getNAtoms();
	int getNActiveAtoms();
	int getNResidues();
	atom getAtom(int);
	atom getAtomOld(int);
	void getActiveAtoms(atom *a);
	void getActiveAtomsOld(atom *a);

	void connectAtoms();
	vectorOfIntegers getConnectivity(int);

	void assignMasses();

	double getMassOfAtom(int);
		
	void expressXYZinMF();
	
	void setAtomPosition(int,vectorOfDoubles);
	vectorOfDoubles getAtomPosition(int);

	void setCell(int,int,int,int);
	
	frameRT getMF();
	
	void calculateCenterOfMass();
	vectorOfDoubles getCenterOfMass();

	void calculateInertiaTensor();
	void calculateInertiaTensorInCM();
	matrix3D getInertiaTensor();

	vectorOfDoubles applyRTAtoms(std::string,frameRT);
	int getNumberOfAtoms(std::string);
	vectorOfIntegers getIndexesOfAtoms(std::string);

	void scaleAtomsPositions(double);

	// Z-Matrix related stuff
	
	void setMainDihedralAngle(int, int, int, int);
	void buildZMatrix(void);
	void outputZMatrix(std::string);
	//void loadZMatrix(std::string); // TODO
	void changeZMatrix(int, int, double, int, double, int, double);
	void changeZMatrix(int, zmatline);
	void buildXYZfromZMatrix(void);
	vectorOfDoubles getFirstDerivative(int, int, int);
	vectorOfDoubles getSecondDerivative(int, int, int, int, int);
	void setZMatrix(zmatline* z); // WARNING: this should be used only for debugging purposes

	// Methods related to calculation of position and derviatives of atoms expressed in MF
	void calculateDerivativesOfCM(void);
	vectorOfDoubles getAtomPositionInMF(int atID); // atID is 1-based
	vectorOfDoubles getFirstDerivativeInMF(int a1, int q, int a2);
	vectorOfDoubles getSecondDerivativeInMF(int a1, int q2, int a2, int q3, int a3);

	// Check if refence atoms have been set
	int areRefAtomsSet(void);

	// Check if Z-Mat has been built
	int hasZMatrix(void);

	// Custom roto-translation
	void setCustomOrientation(double al, double be, double ga);
	void setCustomOrientation(int a1, int a2, int a3);
	void setCustomOrientation(matrix3D m);
	void setCustomCenter(vectorOfDoubles T);
	void setCustomCenter(int a1);
	vectorOfDoubles getAtomPositionInCustomFrame(int atID);
	vectorOfDoubles getFirstDerivativeInCustomFrame(int a1, int q, int a2);
	vectorOfDoubles getSecondDerivativeInCustomFrame(int a1, int q2, int a2, int q3, int a3);

private:
	bool fileIsSet;
	bool moleculeIsLoaded;
	int nAtoms;
	int nActiveAtoms;
	int nResidues;
	vectorOfIntegers CarbonAtoms;
	vectorOfIntegers NitrogenAtoms;
	vectorOfIntegers OxygenAtoms;
	vectorOfIntegers OtherAtoms;
	vectorOfDoubles centerOfMass; // position of CM
	vectorOfDoubles dRCM; // first derivatives of CM
	vectorOfDoubles d2RCM; // second derivatives of CM
	vectorOfDoubles orientation;
	matrix3D inertiaTensor;
	matrix3D rotationToMF; // matrix to rotate from BF to CF (== MF)
	std::string fileName;
	vectorOfAtoms atoms;
	config *conf;
	frameRT MF;
	
	int readXYZ();
	int saveXYZ();
	int readPDB();
	int savePDB();
	double getRadiusOfAtom(int,vectorOfStrings,vectorOfDoubles);
	jacobi diag;
	
	// Z-Matrix related stuff
	
	int isZMatrixBuilt;

	int dA1, dA2, dA3, dA4;
	void fragmentMolecule(void);
	void buildAmat(void);
	void permute(void);
	double getBondLength(int, int);
	double getBondAngle(int, int, int);
	double getDihedralAngle(int, int, int, int);
	int minBond(int, vectorOfIntegers);

	// Fragments data

	int referenceAtomsAreSet;

	vectorOfIntegers atomsID; // 0-based ID's of atoms
	vectorOfIntegers atomsActiveID; // 0-based ID's of active atoms
	int n_dihedrals, nfragments;                            
	int *nterm;                                                                 
	int *ordering;                                                  
	frag *fragments;                                                            

	// Support routines
	inline int qidx(int qid);

	// Custom roto-translation
	matrix3D customRot;
	vectorOfDoubles customCenter;

};

#endif
