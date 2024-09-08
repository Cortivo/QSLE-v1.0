#include "molecule.h"

/****************/
/* Constructors */
/****************/

molecule::molecule()
{
	nAtoms = 0;
	centerOfMass = vectorOfDoubles(3);
	orientation = vectorOfDoubles(3);
	moleculeIsLoaded = false;
	fileIsSet = false;
	fileName = "";
	atoms.clear();

	// Some checks
	referenceAtomsAreSet = 0;
	isZMatrixBuilt = 0;

	// other inits
	customRot.xx = 1.0; customRot.xy = 0.0; customRot.xz = 0.0;
	customRot.yx = 0.0; customRot.yy = 1.0; customRot.yz = 0.0;
	customRot.zx = 0.0; customRot.zy = 0.0; customRot.zz = 1.0;
	customCenter = vectorOfDoubles(3, 0.0);

	return;
}

molecule::molecule(config *c)
{
	nAtoms = 0;
	centerOfMass = vectorOfDoubles(3);
	orientation = vectorOfDoubles(3);
	moleculeIsLoaded = false;
	fileIsSet = false;
	fileName = "";
	atoms.clear();
	conf = c;

	// Some checks
	referenceAtomsAreSet = 0;
	isZMatrixBuilt = 0;

	// other inits
	customRot.xx = 1.0; customRot.xy = 0.0; customRot.xz = 0.0;
	customRot.yx = 0.0; customRot.yy = 1.0; customRot.yz = 0.0;
	customRot.zx = 0.0; customRot.zy = 0.0; customRot.zz = 1.0;
	customCenter = vectorOfDoubles(3, 0.0);

	return;
}

/*****************/
/* Deconstructor */
/*****************/
molecule::~molecule()
{
	atoms.clear();
	return;
}

/********************************************************/
/* Set the file name for I/O operation of this molecule */
/********************************************************/
void molecule::setIOFileName(std::string f)
{
	fileName = f;
	fileIsSet = true;
	#ifdef DEBUG
	std::cout << "INFO >> PDB file name set to : " << fileName << std::endl;
	#endif
	return;
}

/**************************************/
/* Return the fileName, if it was set */
/**************************************/
std::string molecule::getFileName()
{
	if (fileIsSet)
		return fileName;
	else
	{
		std::cout << "ERROR >> molecule::getFileName() : I/O file name did not set yet. Returning NULL string" << std::endl;
		return NULL;
	}
}

/***********************************/
/* Read the molecule from fileName */
/***********************************/
int molecule::loadMolecule()
{
	std::fstream iFile;
	int operationCode = 0; // 0 = OK, 1 = No file name set , 2 = cannot open file in read-mode
	if (!fileIsSet)
	{
		std::cout << "ERROR >> PDB file name is not set. Returning error code 1" << std::endl;
		operationCode = 1;
	}
	else
	{
		int dotPos = fileName.find_last_of(".",fileName.length()) + 1;
		std::string extension = fileName.substr(dotPos,fileName.length());
		#ifdef DEBUG
		std::cout << "INFO >> molecule::loadMolecule() : file extention is : " << extension << std::endl;
		#endif
		if (!extension.compare("pdb") | !extension.compare("PDB"))
			readPDB();
		else
		{
			std::cout << "ERROR >> extension " << extension << " not implemented yet; please convert to pdb" << std::endl;
			exit(1);
		}
		#ifdef DEBUG
		std::cout << "INFO >> successfully loaded molecule from file " << fileName << std::endl;
		#endif
		moleculeIsLoaded = true;
		operationCode = 0;
	}
	
	return operationCode;
}

/*********************************/
/* Save the molecule to fileName */
/*********************************/
int molecule::dumpMolecule()
{
	std::fstream oFile;
	int operationCode = 0; // 0 = OK, 1 = No file name set , 2 = cannot open file in write-mode
	if (!fileIsSet)
	{
		std::cout << "ERROR >> PDB file name is not set. Returning error code 1" << std::endl;
		operationCode = 1;
	}
	else
	{
		int dotPos = fileName.find_last_of(".",fileName.length()) + 1;
		std::string extension = fileName.substr(dotPos,fileName.length());
		#ifdef DEBUG
		std::cout << "INFO >> molecule::dumpMolecule() : file extention is : " << extension << std::endl;
		if (!extension.compare("xyz") | !extension.compare("XYZ"))
			saveXYZ();
		else if (!extension.compare("pdb") | !extension.compare("PDB"))
			savePDB();
		else
		{
			std::cout << "ERROR >> extension " << extension << " not recognized; please check syntax" << std::endl;
			exit(1);
		}
		std::cout << "INFO >> molecule::dumpMolecule() : successfully saved molecule to file " << fileName << std::endl;
		#endif
		operationCode = 0;
	}
	
	return operationCode;
}

/******************************/
/* Return the number of atoms */
/******************************/
int molecule::getNAtoms()
{
	if (!moleculeIsLoaded)
		std::cout << "WARNING >> molecule::getNAtoms() : number of atoms is 0 because no molecule has been loaded yet via loadMolecule() method" << std::endl;
	return nAtoms;
}

int molecule::getNActiveAtoms()
{
	if (!moleculeIsLoaded)
		std::cout << "WARNING >> molecule::getNAtoms() : number of atoms is 0 because no molecule has been loaded yet via loadMolecule() method" << std::endl;
	return nActiveAtoms;
}

/*********************************/
/* Return the number of residues */
/*********************************/
int molecule::getNResidues()
{
	if (!moleculeIsLoaded)
		std::cout << "WARNING >> molecule::getNAtoms() : number of atoms is 0 because no molecule has been loaded yet via loadMolecule() method" << std::endl;
	return nResidues;
}

/********************/
/* Return i-th atom */
/********************/
// i is 1-based index in the NEW numeration
atom molecule::getAtom(int i)
{
	if (!moleculeIsLoaded)
	{
		std::cout << "ERROR >> molecule::getAtom(int) : no molecule has been loaded yet via loadMolecule() method" << std::endl;
		exit(1);
	}

	return atoms[atomsID[i-1]];
}

// i is 1-based index in the OLD numeration
atom molecule::getAtomOld(int i)
{
	if (!moleculeIsLoaded)
	{
		std::cout << "ERROR >> molecule::getAtom(int) : no molecule has been loaded yet via loadMolecule() method" << std::endl;
		exit(1);
	}

	return atoms[i-1];
}

/***********************************/
/* Return the list of active atoms */
/***********************************/
void molecule::getActiveAtoms(atom *a)
{
	if (!moleculeIsLoaded)
	{
		std::cout << "ERROR >> molecule::getAtom(int) : no molecule has been loaded yet via loadMolecule() method" << std::endl;
		exit(1);
	}
	int ia = 0;
	for (int ja = 0; ja < nAtoms; ++ja)
	{
		if (atoms[atomsID[ja]].isActive())
		{
			a[ia] = atoms[atomsID[ja]];
			ia++;
		}
	}
	return;
}

void molecule::getActiveAtomsOld(atom *a)
{
	if (!moleculeIsLoaded)
	{
		std::cout << "ERROR >> molecule::getAtom(int) : no molecule has been loaded yet via loadMolecule() method" << std::endl;
		exit(1);
	}
	int ia = 0;
	for (int ja = 0; ja < nAtoms; ++ja)
	{
		if (atoms[ja].isActive())
		{
			a[ia] = atoms[ja];
			ia++;
		}
	}
	return;
}

/******************************/
/* Echo molecule informations */
/******************************/
std::string molecule::toString(bool printID, bool printAtomTypes, bool printCoordinates, bool printConnectivity, bool printResidueName, bool printResidueNumber)
{
	std::ostringstream molInfo;
	vectorOfIntegers b;
	molInfo << "DATA OF MOLECULE " << fileName << std::endl;
	for (int i = 0; i < nAtoms; i++)
	{
		if (printID)
			molInfo << (i+1) << ")\t";
		if (printAtomTypes)
			molInfo << atoms.at(i).getAtomType() << "\t";
		if (printCoordinates)
			molInfo << std::fixed << std::setprecision(3) << atoms.at(i).x() << "\t" << atoms.at(i).y() << "\t" << atoms.at(i).z() << "\t";
		if (printConnectivity)
		{
			b = atoms.at(i).getBonds();
			for (unsigned int j = 0; j < b.size(); j++)
				molInfo << b.at(j) << "\t";
		}
		if (printResidueName)
			molInfo << atoms.at(i).getResidueName() << "\t";
		if (printResidueNumber)
			molInfo << atoms.at(i).getResidueNumber() << "\t";
		molInfo << std::endl;
	}
	return molInfo.str();
}

/**********************/
/* XYZ I/O operations */
/**********************/

int molecule::readXYZ()
{
	int operationCode = 0;
	
	/* 1. Convert molecule to PDB using SeStO */
	
	std::fstream iFile;
	iFile.open(fileName.c_str(),std::ios::in);
	if (!iFile.is_open())
	{
		std::cout << "ERROR >> molecule::readXYZ() : cannot open/read file " << fileName << std::endl;
		exit(1);
	}
	iFile.close();

	std::cout << "INFO >> molecule::readXYZ() : converting XYZ file " << fileName << " to PDB using SeStO  - remember to inspect the sesto_load.log file" << std::endl;

	int dotPos = fileName.find_last_of(".",fileName.length()) + 1;
	int slashPos = fileName.find_last_of("/",fileName.length()) + 1;
	std::string command = conf->getSesto() + " " + fileName.substr(0,slashPos) + " " + fileName.substr(slashPos,dotPos-slashPos-1) + \
					 " " + conf->getVdwFile() +" " + conf->getAminoAcidsFile() + " > ./sesto_load.log";
	//system(command.c_str());
	
	/* 2. Read the PDB file */
	
	std::string extension = fileName.substr(dotPos,fileName.length());
	fileName.replace(dotPos, 3, "pdb");
	readPDB();

	return operationCode;	
}

int molecule::saveXYZ()
{
	std::fstream oFile;
	oFile.open(fileName.c_str(),std::ios::out);
	if (!oFile.is_open())
	{
		std::cout << "ERROR >> molecule::saveXYZ() : cannot create/write file " << fileName << std::endl;
		exit(1);
	}
	oFile << nAtoms << std::endl;
	oFile << fileName << std::endl;
	for (int i = 0; i < nAtoms; i++)
		oFile << atoms[atomsID[i]].getAtomType() << "\t" << atoms[atomsID[i]].x() << "\t" << atoms[atomsID[i]].y() << "\t" << atoms[atomsID[i]].z() << std::endl;
	oFile.close();
	return 0;
}

/**********************/
/* PDB I/O operations */
/**********************/

int molecule::readPDB()
{
	int operationCode = 0;

	int stop;

	CarbonAtoms.clear();
	NitrogenAtoms.clear();
	OxygenAtoms.clear();
	OtherAtoms.clear();
	
	#ifdef DEBUG
	std::cout << "INFO >> reading PDB file " << fileName << std::endl;
	#endif
	
	std::fstream iFile;
	iFile.open(fileName.c_str(),std::ios::in);
	if (!iFile.is_open())
	{
		std::cout << "ERROR >> cannot read/open PDB file " << fileName << std::endl;
		exit(1);
	}
	char fileLine[2048];
	std::string fileStr;
	std::string tmpStr;
	int spacePos;
	int pos1, pos2;
	
	int residueNumber;
	int isActive, activeID = 1;
	double x, y, z, tmpDouble;
	std::string atomType;
	std::string residueName;
	std::string element;

	nAtoms = 0;
	nActiveAtoms = 0;
	nResidues = 0;
	atoms.clear();
	atom tmpAtom;
	
	atomsID.clear();
	atomsActiveID.clear();

	while (iFile.getline(fileLine,2048))
	{
		fileStr.assign(fileLine);
		spacePos = fileStr.find_first_of(" ",0);
		if (spacePos >= 0)
		{
			/*************************************/
			/* Check for the HETATM or ATOM card */
			/*************************************/
			if (!fileStr.substr(0,spacePos).compare("HETATM") | !fileStr.substr(0,spacePos).compare("ATOM"))
			{
				
				/**************************/
				/* Read the PDB atom type */
				/**************************/
				
				spacePos = fileStr.find(' ',0);
				pos1 = fileStr.find_first_not_of(' ',spacePos+1);
				pos2 = fileStr.find(' ',pos1+1); // found atom ID
				pos1 = fileStr.find_first_not_of(' ',pos2+1);
				pos2 = fileStr.find(' ',pos1+1);
				atomType = fileStr.substr(pos1,pos2-pos1+1);
				
				/*************************/
				/* Read the residue name */
				/*************************/
				
				pos1 = fileStr.find_first_not_of(' ',pos2+1);
				pos2 = fileStr.find(' ',pos1+1);
				residueName = fileStr.substr(pos1,3);

				/***************************/
				/* Read the residue number */
				/***************************/
				
				pos1 = fileStr.find_first_not_of(' ',pos2+1);
				pos2 = fileStr.find(' ',pos1+1);
				tmpStr = fileStr.substr(pos1,pos2-pos1+1);
				if (!isdigit(tmpStr[0]))
				{
					pos1 = fileStr.find_first_not_of(' ',pos2+1);
					pos2 = fileStr.find(' ',pos1+1);
					tmpStr = fileStr.substr(pos1,pos2-pos1+1);					
				}
				sscanf(tmpStr.c_str(),"%d",&residueNumber);
				
				/********************/
				/* Read coordinates */
				/********************/
				
				pos1 = fileStr.find_first_not_of(' ',pos2+1);
				pos2 = fileStr.find(' ',pos1+1);
				tmpStr = fileStr.substr(pos1,pos2-pos1+1);
				sscanf (tmpStr.c_str(),"%lf",&x);
				pos1 = fileStr.find_first_not_of(' ',pos2+1);
				pos2 = fileStr.find(' ',pos1+1);
				tmpStr = fileStr.substr(pos1,pos2-pos1+1);
				sscanf (tmpStr.c_str(),"%lf",&y);
				pos1 = fileStr.find_first_not_of(' ',pos2+1);
				pos2 = fileStr.find(' ',pos1+1);
				tmpStr = fileStr.substr(pos1,pos2-pos1+1);
				sscanf (tmpStr.c_str(),"%lf",&z);
				
				/******************/
				/* Read occupancy */
				/******************/

				pos1 = fileStr.find_first_not_of(' ',pos2+1);
				pos2 = fileStr.find(' ',pos1+1);
				tmpStr = fileStr.substr(pos1,pos2-pos1+1);
				sscanf (tmpStr.c_str(),"%lf",&tmpDouble);
				isActive = (int)tmpDouble;
				nActiveAtoms += isActive;

				/*************/
				/* Read Atom */
				/*************/

				stop = 0;
				if (fileStr.length() < 78)
					stop = 1;
				else
				{
					element = fileStr.substr(76,77);
					element.erase(remove(element.begin(), element.end(), ' '), element.end());
					if (element.empty())
						stop = 1;
				}

				if (stop)
				{
					std::cout << std::endl << std::endl;
					std::cout << "ERROR - No element found at columns 77-78 of atom " << nAtoms+1 << " in the PDB file" << std::endl;
					std::cout << std::endl << std::endl;
					exit(1);
				}

				/**********************/
				/* Update atoms array */
				/**********************/
				
				tmpAtom.x(x);
				tmpAtom.y(y);
				tmpAtom.z(z);
				tmpAtom.isActive(isActive);
				///////////tmpAtom.setAtomTypePDB(atomType); -- TO CANCEL
				tmpAtom.setAtomType(element);
				tmpAtom.setResidueName(residueName);
				tmpAtom.setResidueNumber(residueNumber);
				tmpAtom.clearBonds();
				atoms.push_back(tmpAtom);
				atoms[nAtoms].id = atoms[nAtoms].id0 = nAtoms + 1;
				if (isActive)
				{
					atoms[nAtoms].activeID = atoms[nAtoms].activeID0 = activeID;
					atomsActiveID.push_back(nActiveAtoms);
					activeID++;
				}
				atomsID.push_back(nAtoms);
				atoms[nAtoms].savePosition();
				if (nAtoms == 0)
					nResidues = 1;
				else if (residueNumber != atoms[nAtoms - 1].getResidueNumber())
					nResidues++;
				
				if (!atoms[nAtoms].getAtomType().compare("C") | !atoms[nAtoms].getAtomType().compare("c"))
					CarbonAtoms.push_back(nAtoms);
				else if (!atoms[nAtoms].getAtomType().compare("N") | !atoms[nAtoms].getAtomType().compare("n"))
					NitrogenAtoms.push_back(nAtoms);
				else if (!atoms[nAtoms].getAtomType().compare("O") | !atoms[nAtoms].getAtomType().compare("o"))
					OxygenAtoms.push_back(nAtoms);
				else
					OtherAtoms.push_back(nAtoms);

				nAtoms++;
				
			}
		}
	}
	#ifdef DEBUG
	std::cout << "INFO >> molecule::readPDB() : Found " << nActiveAtoms << " active atoms in the PDB. Thus, the number of free internal degrees of freedom will be " << 3 * nActiveAtoms - 6 << ". " << std::endl;
	std::cout << "INFO >> PDB file succesfully read" << std::endl;
	#endif

	/***************************/
	/* Estabilish connectivity */
	/***************************/

	connectAtoms();

	/*****************/
	/* Assign masses */
	/*****************/
	
	assignMasses();

	return operationCode;		
}

int molecule::savePDB()
{
	int operationCode = 0;
	
	// 1. Write file in XYZ format 
	int dotPos = fileName.find_last_of(".",fileName.length()) + 1;
	fileName.replace(dotPos, 3, "xyz");
	saveXYZ();
	
	// 2. Convert file to PDB using SeStO

	std::cout << "INFO >> molecule::savePDB() : converting XYZ file " << fileName << " to PDB using SeStO - remember to inspect the sesto_save.log file" << std::endl;
	int slashPos = fileName.find_last_of("/",fileName.length()) + 1;
	std::string command = conf->getSesto() + " " + fileName.substr(0,slashPos) + " " + fileName.substr(slashPos,dotPos-slashPos-1) + \
					 " " + conf->getVdwFile() +" " + conf->getAminoAcidsFile() + " > ./sesto_save.log";
	//system(command.c_str());	
	
	return operationCode;
}

/*********************************/
/* Find connectivity among atoms */
/*********************************/
void molecule::connectAtoms()
{
	#ifdef DEBUG
	std::cout << "INFO >> determining connectivity between atoms" << std::endl;
	#endif	
	/*******************************/
	/* Read the radii file vdw.dat */
	/*******************************/
	
	std::fstream radiiFile;
	radiiFile.open(conf->getVdwFile().c_str(),std::ios::in);
	char fileLine[2048];
	char tmpChr[4];
	std::string fileStr;
	std::string tmpStr;
	double tmpDouble;
	
	vectorOfStrings types;
	vectorOfDoubles radii;
	
	while (radiiFile.getline(fileLine,2048))
	{
		fileStr.assign(fileLine);
		if (fileStr.substr(0,1).compare("#") & fileStr.substr(0,1).compare(" ")) // lines starting with # or a white space are considered as comments
		{
			sscanf(fileStr.c_str(),"%s %lf",tmpChr,&tmpDouble);
			tmpStr.assign(tmpChr);
			types.push_back(tmpStr);
			radii.push_back(tmpDouble);
		}
	}
	int nRadii = radii.size();
	

	#ifdef DEBUG
	std::cout << std::endl << "Atom\tRadius / A" <<std::endl << "------------------" << std::endl;
	for (int i = 0; i < nRadii; i++) std::cout << types.at(i) << "\t" << radii.at(i) << std::endl;
	std::cout << "------------------" << std::endl << std::endl;
	#endif
	
	/*****************/
	/* Connect atoms */
	/*****************/
	
	double xi, yi, zi;
	double xj, yj, zj;
	double rij, r_vdw_ij;
	
	for (int i = 0; i < nAtoms; i++)
	{
		xi = atoms[i].x(); yi = atoms[i].y(); zi = atoms[i].z();
		for (int j = 0; j < nAtoms; j++)
		{
			if (j != i)
			{
				xj = atoms[j].x(); yj = atoms[j].y(); zj = atoms[j].z();
				r_vdw_ij = getRadiusOfAtom(i,types,radii) + getRadiusOfAtom(j,types,radii) + SLOP_FACTOR;
				r_vdw_ij *= r_vdw_ij;
				rij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) + (zj-zi)*(zj-zi);
				if (rij <= r_vdw_ij)
					atoms.at(i).addBond(j+1);
			}
		}
		std::string aname = atoms[i].getAtomType();
		if (!aname.compare("H") && atoms[i].getNBonds() > 1)
		{
			std::cout << std::endl << "ERROR: hydrogen atom ID " << i+1 << " seems to be linked to more than 1 atom" << std::endl;
			exit(1);
		}
	}
	#ifdef DEBUG	
	std::cout << "INFO >> connectivity done" << std::endl;
	#endif	
	return;
}

/*****************************************/
/* Return the connectivity of given atom */
/*****************************************/
// atom id 1-based, permuted order
vectorOfIntegers molecule::getConnectivity(int aid)
{
	vectorOfIntegers c0, c;
	int a = aid - 1;

	c0 = atoms[atomsID[a]].getBonds();

	for (int i = 0; i < c0.size(); ++i)
		c.push_back(atoms[c0[i]-1].id);

	return c;
}

/******************************************************/
/* Return the raidus of a-th atom. Index a is 0-based */
/******************************************************/
double molecule::getRadiusOfAtom(int a, vectorOfStrings types, vectorOfDoubles radii)
{
	for (unsigned int i = 0; i < radii.size(); i++)
	{
		if (!atoms.at(a).getAtomType().compare(types.at(i)))
			return radii.at(i);
	}
	std::cout << "ERROR >> molecule::getRadiusOfAtom(int, vectorOfStrings, vectorOfDoubles) : cannot find radius for element " << atoms.at(a).getAtomType() << std::endl;
	exit(1);
	return 0.0;
}

/**************************/
/* Assign masses to atoms */
/**************************/
void molecule::assignMasses()
{
	std::fstream massFile;
	massFile.open(conf->getMassFile().c_str(),std::ios::in);
	char fileLine[2048];
	char tmpChr[4];
	std::string fileStr;
	std::string tmpStr;
	double tmpDouble;
	
	vectorOfStrings types;
	vectorOfDoubles mass;
	
	while (massFile.getline(fileLine,2048))
	{
		fileStr.assign(fileLine);
		if (fileStr.substr(0,1).compare("#") & fileStr.substr(0,1).compare(" ")) // lines starting with # or a white space are considered as comments
		{
			sscanf(fileStr.c_str(),"%s %lf",tmpChr,&tmpDouble);
			tmpStr.assign(tmpChr);
			types.push_back(tmpStr);
			mass.push_back(tmpDouble);
		}
	}
	int nMass = mass.size();
	#ifdef DEBUG
	std::cout << "INFO >> molecule::assignMasses() : Using the following set of masses" << std::endl;
	std::cout << std::endl << "Atom\tMass / uma" <<std::endl << "------------------" << std::endl;
	for (int i = 0; i < nMass; i++) std::cout << types.at(i) << "\t" << mass.at(i) << std::endl;
	std::cout << "------------------" << std::endl << std::endl;
	#endif
	for (int i = 0; i < nAtoms; i++)
		atoms[i].setMass(types,mass);
}

/****************************************************/
/* Return the mass of a-th atom. Index a is 0-based */
/****************************************************/
double molecule::getMassOfAtom(int a)
{
	return atoms[a].getMass();
}

/*******************************************/
/* Express the atoms coordinates in the MF */
/*******************************************/
void molecule::expressXYZinMF()
{	
	vectorOfDoubles zero3 = vectorOfDoubles(3,0.0);
	matrix3D identity3 = matrix3D(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
	

	// Center of mass
	calculateCenterOfMass();

	// translate coordinates
	MF.setTranslation(centerOfMass);
	MF.setRotation(identity3);
	for (int i = 0; i < nAtoms; i++)
		atoms[i].setPosition(MF.transform(atoms[i].getPosition()));

	// calculate and diagonalize inertia tensor
	calculateInertiaTensor();
	diag.setMaxStep(500);
	diag.setMatrix(inertiaTensor);
	diag.diagonalize();
	diag.reorder();
	/**diag.debug();**/

	// Rotate molecule
	matrix3D Emat = diag.getEigenVectors3D();
	Emat.transpose();
	MF.setTranslation(zero3);
	MF.setRotation(Emat);
	for (int i = 0; i < nAtoms; i++)
		atoms[i].setPosition(MF.transform(atoms[i].getPosition()));

	// update center of mass and inertia tensor
	calculateCenterOfMass();
	calculateInertiaTensor();

	return;
}

/********************************/
/* Calculate the center of mass */
/********************************/
void molecule::calculateCenterOfMass()
{
	centerOfMass = vectorOfDoubles(3,0.0);
	double xi,yi,zi;
	double mass, totMass = 0.0;
	double cmx = 0.0, cmy = 0.0, cmz = 0.0;
	for (int i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].getMass();
		xi = mass*atoms[i].x();
		yi = mass*atoms[i].y();
		zi = mass*atoms[i].z();
		cmx += xi;
		cmy += yi;
		cmz += zi;
		totMass += mass;
	}
	cmx /= totMass;
	cmy /= totMass;
	cmz /= totMass;
	centerOfMass[0] = cmx;
	centerOfMass[1] = cmy;
	centerOfMass[2] = cmz;

	return;
}
vectorOfDoubles molecule::getCenterOfMass()
{
	return centerOfMass;
}

/******************/
/* Inertia tensor */
/******************/
void molecule::calculateInertiaTensor()
{
	double xi,yi,zi;
	double mass;
	double Ixx = 0.0, Ixy = 0.0, Ixz = 0.0;
	double Iyx = 0.0, Iyy = 0.0, Iyz = 0.0;
	double Izx = 0.0, Izy = 0.0, Izz = 0.0;
	for (int i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].getMass();
		xi = atoms[i].x();
		yi = atoms[i].y();;
		zi = atoms[i].z();

		Ixx += mass * (yi * yi + zi * zi);
		Ixy -= mass * xi * yi;
		Ixz -= mass * xi * zi;
		
		Iyy += mass * (zi * zi + xi * xi);
		Iyz -= mass * yi * zi;

		Izz += mass * (xi * xi + yi * yi);

	}

	Iyx = Ixy;
	Izx = Ixz;
	Izy = Iyz;
	inertiaTensor = matrix3D(Ixx,Ixy,Ixz,Iyx,Iyy,Iyz,Izx,Izy,Izz);
	return;
}

matrix3D molecule::getInertiaTensor()
{
	return inertiaTensor;
}

void molecule::calculateInertiaTensorInCM()
{
	double xi,yi,zi;
	double mass;
	double Ixx = 0.0, Ixy = 0.0, Ixz = 0.0;
	double Iyx = 0.0, Iyy = 0.0, Iyz = 0.0;
	double Izx = 0.0, Izy = 0.0, Izz = 0.0;

	calculateCenterOfMass();
	vectorOfDoubles rcm = getCenterOfMass();

	for (int i = 0; i < nAtoms; i++)
	{
		mass = atoms[i].getMass();
		xi = atoms[i].x() - rcm[0];
		yi = atoms[i].y() - rcm[1];
		zi = atoms[i].z() - rcm[2];

		Ixx += mass * (yi * yi + zi * zi);
		Ixy -= mass * xi * yi;
		Ixz -= mass * xi * zi;
		
		Iyy += mass * (zi * zi + xi * xi);
		Iyz -= mass * yi * zi;

		Izz += mass * (xi * xi + yi * yi);

	}

	Iyx = Ixy;
	Izx = Ixz;
	Izy = Iyz;
	inertiaTensor = matrix3D(Ixx,Ixy,Ixz,Iyx,Iyy,Iyz,Izx,Izy,Izz);
	return;
}

/********************/
/* Set cell numbers */
/********************/

void molecule::setCell(int i, int g1, int g2, int g3)
{
	atoms[i-1].setCellNumbers(g1,g2,g3);
	return;
}

/*****************/
/* Return the MF */
/*****************/

frameRT molecule::getMF()
{
	return MF;
}

/***************************************/
/* Modify position vector of i-th atom */
/***************************************/
// i is 1-based counter
void molecule::setAtomPosition(int i, vectorOfDoubles r)
{
	atoms[atomsID[i-1]].setPosition(r);
	return;
}
// a is 1-based counter
vectorOfDoubles molecule::getAtomPosition(int a)
{
	return atoms[atomsID[a-1]].getPosition();
}
/***************************************************************************/
/* Returns rotated (NO PERMANENT MODIFICATION) C, N or O atoms coordinates */
/***************************************************************************/

vectorOfDoubles molecule::applyRTAtoms(std::string AT, frameRT RT)
{
	if (!AT.compare("C") | ! AT.compare("c"))
	{
		int nC = CarbonAtoms.size();
		vectorOfDoubles newCarbonCoordinates = vectorOfDoubles(3*nC,0.0);
		vectorOfDoubles Ccoord = vectorOfDoubles(3,0.0);
		for (int i = 0; i < nC; i++)
		{
			Ccoord = atoms[CarbonAtoms[i]].getPosition();
			Ccoord = RT.transform(Ccoord);
			newCarbonCoordinates[i*3+0] = Ccoord[0];
			newCarbonCoordinates[i*3+1] = Ccoord[1];
			newCarbonCoordinates[i*3+2] = Ccoord[2];
		}
		return newCarbonCoordinates;
	}
	else if (!AT.compare("N") | !AT.compare("n"))
	{
		int nN = NitrogenAtoms.size();
		vectorOfDoubles newNitrogenCoordinates = vectorOfDoubles(3*nN,0.0);
		vectorOfDoubles Ncoord = vectorOfDoubles(3,0.0);
		for (int i = 0; i < nN; i++)
		{
			Ncoord = atoms[NitrogenAtoms[i]].getPosition();
			Ncoord = RT.transform(Ncoord);
			newNitrogenCoordinates[i*3+0] = Ncoord[0];
			newNitrogenCoordinates[i*3+1] = Ncoord[1];
			newNitrogenCoordinates[i*3+2] = Ncoord[2];
		}
		return newNitrogenCoordinates;
	}
	else if (!AT.compare("O") | !AT.compare("o"))
	{
		int nO = OxygenAtoms.size();
		vectorOfDoubles newOxygenCoordinates = vectorOfDoubles(3*nO,0.0);
		vectorOfDoubles Ocoord = vectorOfDoubles(3,0.0);
		for (int i = 0; i < nO; i++)
		{
			Ocoord = atoms[OxygenAtoms[i]].getPosition();
			Ocoord = RT.transform(Ocoord);
			newOxygenCoordinates[i*3+0] = Ocoord[0];
			newOxygenCoordinates[i*3+1] = Ocoord[1];
			newOxygenCoordinates[i*3+2] = Ocoord[2];
		}
		return newOxygenCoordinates;
	}
	else if (!AT.compare("other"))
	{
		int nE = OtherAtoms.size();
		vectorOfDoubles newOtherCoordinates = vectorOfDoubles(3*nE,0.0);
		vectorOfDoubles Ecoord = vectorOfDoubles(3,0.0);
		for (int i = 0; i < nE; i++)
		{
			Ecoord = atoms[OtherAtoms[i]].getPosition();
			Ecoord = RT.transform(Ecoord);
			newOtherCoordinates[i*3+0] = Ecoord[0];
			newOtherCoordinates[i*3+1] = Ecoord[1];
			newOtherCoordinates[i*3+2] = Ecoord[2];
		}
		return newOtherCoordinates;
	}
	else
	{
		std::cout << "ERROR >> molecule::applyRotationToAtoms(std::string) : atom type " << AT << " not recognized. Please choose among C, N or O" << std::endl;
		exit(1);
	}
}
/************************************************/
/* Return the number of atoms of specified type */
/************************************************/
int molecule::getNumberOfAtoms(std::string AT)
{
	     if (!AT.compare("C") | ! AT.compare("c")) return CarbonAtoms.size();
	else if (!AT.compare("N") | ! AT.compare("n")) return NitrogenAtoms.size();
	else if (!AT.compare("O") | ! AT.compare("o")) return OxygenAtoms.size();
	else if (!AT.compare("other"))                 return OtherAtoms.size();
        else
	{
		std::cout << "ERROR >> molecule::getNumberOfAtoms(std::string) : atom type " << AT << " not recognized. Please choose among C, N or O" << std::endl;
		exit(1);
	}
}
/**********************************************************/
/* Return the array of indexes of atoms of specified type */
/**********************************************************/
vectorOfIntegers molecule::getIndexesOfAtoms(std::string AT)
{
	     if (!AT.compare("C") | ! AT.compare("c")) return CarbonAtoms;
	else if (!AT.compare("N") | ! AT.compare("n")) return NitrogenAtoms;
	else if (!AT.compare("O") | ! AT.compare("o")) return OxygenAtoms;
	else if (!AT.compare("other")) return OtherAtoms;
        else
	{
		std::cout << "ERROR >> molecule::getIndexesOfAtoms(std::string) : atom type " << AT << " not recognized. Please choose among C, N or O" << std::endl;
		exit(1);
	}
}

/**************************/
/* Z-Matrix related stuff */
/**************************/
void molecule::setMainDihedralAngle(int i1, int i2, int i3, int i4)
{
	dA1 = i1;
	dA2 = i2;
	dA3 = i3;
	dA4 = i4;
	n_dihedrals = 1;
	nfragments  = 2;
	fragmentMolecule();
	referenceAtomsAreSet = 1;
	return;
}


void molecule::buildZMatrix(void)
{
	zmatline z;

	// Reorder atoms to build Z-Matrix

	permute();
	#ifdef DEBUG
	std::cout << "Applied permutation:"<< std::endl << std::endl;
	std::cout << "old ID  | new ID" << std::endl;
	std::cout << "-----------------" << std::endl;
	for (int i = 0; i < nAtoms; ++i)
		std::cout << atoms[i].id0 << "\t  " << atoms[i].id <<  std::endl;

	std::cout << std::endl;
	
	std::cout << "ID of MODEL (active) atoms:"<< std::endl << std::endl;
	std::cout << "old Model ID  | new Model ID" << std::endl;
	std::cout << "----------------------------" << std::endl;
	for (int i = 0; i < nAtoms; ++i)
	{
		if (atoms[i].isActive())
			std::cout << atoms[i].activeID0 << "\t\t  " << atoms[i].activeID <<  std::endl;
	}
	#endif
	// Build Z-Mat entries for atoms 1, 2, 3, and 4

	atoms[atomsID[0]].setZMatrixEntry(0, 0.0, 0, 0.0, 0, 0.0);
	atoms[atomsID[1]].setZMatrixEntry(1, getBondLength(atomsID[0]+1, atomsID[1]+1), 0, 0.0, 0, 0.0);
	atoms[atomsID[2]].setZMatrixEntry(2, getBondLength(atomsID[1]+1, atomsID[2]+1), 1, getBondAngle(atomsID[0]+1, atomsID[1]+1, atomsID[2]+1), 0, 0.0);
	atoms[atomsID[3]].setZMatrixEntry(3, getBondLength(atomsID[2]+1, atomsID[3]+1), 2, getBondAngle(atomsID[1]+1, atomsID[2]+1, atomsID[3]+1), 1, getDihedralAngle(atomsID[0]+1, atomsID[1]+1, atomsID[2]+1, atomsID[3]+1));

	// This override is to set atom 2 as reference, in the XYZ coordinates
	// Here we are not changing the Z-mat entries; just swapping the Bmat
	// between atoms 1 and 2
	atoms[atomsID[0]].overrideBmat(atoms[atomsID[1]].getZMatrixEntry());
	atoms[atomsID[1]].overrideBmat(atoms[atomsID[0]].getZMatrixEntry());

	// Build Z-Mat entries for all the other atoms

	vectorOfIntegers L1(1);
	vectorOfIntegers L2(2);
	vectorOfIntegers L3(3);

	//std::cout << std::endl << std::endl;
	for (int i = 4; i < nAtoms; ++i)
	{
		L1[0] = atomsID[i]+1;
		z.A = minBond(atomsID[i], L1);
		z.A = atoms[z.A-1].id;
		z.d = getBondLength(atomsID[i]+1, atomsID[z.A-1]+1);

		if (z.A == 2)
		{
			z.B = 1;
			z.theta = getBondAngle(atomsID[i]+1, atomsID[z.A-1]+1, atomsID[z.B-1]+1);

			z.C = 3;
			z.phi = getDihedralAngle(atomsID[i]+1, atomsID[z.A-1]+1, atomsID[z.B-1]+1, atomsID[z.C-1]+1);
		}
		else
		{
			L2[0] = atomsID[i]+1;
			L2[1] = atomsID[z.A-1]+1;
			z.B = minBond(atomsID[z.A-1], L2);
			z.B = atoms[z.B-1].id;
			z.theta = getBondAngle(atomsID[i]+1, atomsID[z.A-1]+1, atomsID[z.B-1]+1);

			L3[0] = atomsID[i]+1;
			L3[1] = atomsID[z.A-1]+1;
			L3[2] = atomsID[z.B-1]+1;
			z.C = minBond(atomsID[z.B-1], L3);
			z.C = atoms[z.C-1].id;
			z.phi = getDihedralAngle(atomsID[i]+1, atomsID[z.A-1]+1, atomsID[z.B-1]+1, atomsID[z.C-1]+1);
		}

		atoms[atomsID[i]].setZMatrixEntry(z);
	}

	outputZMatrix("screen");
	outputZMatrix("zmatrix.zmt");

	buildAmat();

	isZMatrixBuilt = 1;

	return;
}

// WARNING: this should be used only for debugging purposes
void molecule::setZMatrix(zmatline* z)
{
	for (int i = 0; i < nAtoms; ++i)
		atoms[atomsID[i]].setZMatrixEntry(z[i]);
	isZMatrixBuilt = 1;
	return;
}

void molecule::outputZMatrix(std::string s)
{
	zmatline z;
	// Standard output
	if (!s.compare("screen"))
	{
		//std::cout << "Z-Matrix" << std::endl;
		for (int i = 0; i < nAtoms; ++i)
		{
			z = atoms[atomsID[i]].getZMatrixEntry();
			//std::cout << i+1 << "\t" << z.A << "\t" << z.d << "\t" << z.B << "\t" << z.theta * RAD2DEG << "\t" << z.C << "\t" << z.phi * RAD2DEG << std::endl;
		}
		//std::cout << std::endl << std::endl;
	}
	// Output on file
	else
	{
		std::fstream f;
		f.open(s.c_str(), std::ios::out);
		for (int i = 0; i < nAtoms; ++i)
		{
			z = atoms[atomsID[i]].getZMatrixEntry();
			f << /*atoms[atomsID[i]].getAtomType()*/ i+1 << "\t" << z.A << "\t" << z.d << "\t" << z.B << "\t" << z.theta * RAD2DEG << "\t" << z.C << "\t" << z.phi * RAD2DEG << std::endl;
		}
		f.close();
	}
	return;
}

// Methods to modify the Z-Matrix
// ID of atoms are 1-based
void molecule::changeZMatrix(int aid, int A, double d, int B, double theta, int C, double phi)
{
	atoms[atomsID[aid-1]].setZMatrixEntry(A, d, B, theta, C, phi);
	buildXYZfromZMatrix();
	return;
}

void molecule::changeZMatrix(int aid, zmatline z)
{
	atoms[atomsID[aid-1]].setZMatrixEntry(z);
	buildXYZfromZMatrix();
	return;
}

// Re-build cartesian coordinates with respect to
// ''internal'' reference placed on atom 2
void molecule::buildXYZfromZMatrix(void)
{
	int N1 = 0;
	for (int i = 0; i < nAtoms; ++i)
		N1 += atoms[i].getNFrag(); // N1 will contain the number of atoms in Fragment 1

	int c1 = 4;// First atom (0-based) in fragment 0
	int c2 = nAtoms - N1 + 2; // First atom (0-based) in fragment 1

	vectorOfDoubles r0(4), r(4);
	matrix4D A;
	matrix4D R(1.0, 0.0, 0.0, 0.0,
	    0.0, -1.0, 0.0, 0.0,
	    0.0, 0.0, -1.0, 0.0,
	    0.0, 0.0, 0.0, 1.0);

	r0[0] = 0.0; r0[1] = 0.0; r0[2] = 0.0; r0[3] = 1.0;

	vectorOfIntegers chain;
	for (int i = 0; i < nAtoms; ++i)
	{
		A = atoms[atomsID[i]].getAmat();
		r = A.multiply(r0);
		if (i >= c1 && i < c2)
		{
			chain = getAtom(i+1).getChain();
			if(atoms[atomsID[i]].getZMatAtom(1) != 2 && chain[chain.size() - 2] == 0)
				r = R.multiply(r); // Apply a 180 deg x-rotation for atoms ''on the left'' of reference atom 1
		}
		atoms[atomsID[i]].x(r[0]);
		atoms[atomsID[i]].y(r[1]);
		atoms[atomsID[i]].z(r[2]);
	}

	return;
}

// Return first derivative of atom a1 (1-based) position
// with respect to Z-Matrix coordinate
// // q = 1 (d), 2 (theta), 3 (phi)
// of atom a2 (1-based)
vectorOfDoubles molecule::getFirstDerivative(int a1, int q, int a2)
{

	if ((a2 == 1) || (a2 == 2 && q != 1) || (a2 == 3 && q == 3))
	{
		std::cout << std::endl << "ERROR in molecule::getFirstDerivative(int a1, int q, int a2) : no internal coordinate q_a2 defined for q = " << q << " and a2 = " << a2 << std::endl << std::endl;
		vectorOfDoubles fd(3);
		fd[0]=0.0;
		fd[1]=0.0;
		fd[2]=0.0;
		return fd;
		//exit(1);
	}

	int N1 = 0;
	for (int i = 0; i < nAtoms; ++i)
		N1 += atoms[i].getNFrag();

	int c1 = 4;
	int c2 = nAtoms - N1 + 2;

	matrix4D R(1.0,  0.0,  0.0,  0.0,
	    	   0.0, -1.0,  0.0,  0.0,
	   	   0.0,  0.0, -1.0,  0.0,
		   0.0,  0.0,  0.0,  1.0);

	int a1id = a1 - 1;
	int a2id = a2 - 1;

	vectorOfDoubles fd(4);
	vectorOfDoubles r0(4);

	matrix4D A(1.0, 0.0, 0.0, 0.0,
		   0.0, 1.0, 0.0, 0.0,
		   0.0, 0.0, 1.0, 0.0,
		   0.0, 0.0, 0.0, 1.0);

	fd[0] = 0.0; fd[1] = 0.0; fd[2] = 0.0; fd[3] = 0.0;
	r0[0] = 0.0; r0[1] = 0.0; r0[2] = 0.0; r0[3] = 1.0;

	switch (a1)
	{
		case 1:
		{
			if (a2 == 2 && q == 1)
			{
				A = atoms[atomsID[0]].getBmatFirstDerivative(1);
				fd = A.multiply(r0);
			}
			break;
		}
		case 2:
		{
			// All derivatives are [0 0 0]
			// being atom 2 the reference atom
			break;
		}
		case 3:
		{
			if (a2 == 3)
			{
				A = atoms[atomsID[2]].getBmatFirstDerivative(q);
				fd = A.multiply(r0);
			}
			break;
		}
		case 4:
		{
			if (a2 == 3)
			{
				A = atoms[atomsID[3]].getBmat();
				A.multiply(atoms[atomsID[2]].getBmatFirstDerivative(q), A);
				fd = A.multiply(r0);
			}
			else if (a2 == 4)
			{
				A = atoms[atomsID[3]].getBmatFirstDerivative(q);
				A.multiply(atoms[atomsID[2]].getBmat(), A);
				fd = A.multiply(r0);
			}
			break;
		}
		default:
		{
			int a3id, count = 0;
			vectorOfIntegers chain = atoms[atomsID[a1id]].getChain();
			intit it = find(chain.begin(), chain.end(), a2id);
			if (it != chain.end())
			{
				a3id = chain[0];
				while (a3id != a2id)
				{
					if (a3id != a2id)
					{
						A.multiply(atoms[atomsID[a3id]].getBmat(), A);
						count++;
						a3id = chain[count];
					}
				}
				A.multiply(atoms[atomsID[a2id]].getBmatFirstDerivative(q), A);
				if (a2id != 2)	// For a2id = 2 is intended atom 3. Thus, the premultiplication by A2 should not be applied (is a translation to atom 1)
				{
					a3id = atoms[atomsID[a2id]].getZMatAtom(1) - 1;
					A.multiply(atoms[atomsID[a3id]].getAmat(), A);
				}
				fd = A.multiply(r0);
			}
		}
	}

	vectorOfIntegers chain;
	if (a1id >= c1 && a1id < c2)
	{
		chain = getAtom(a1id+1).getChain();
		if(atoms[atomsID[a1id]].getZMatAtom(1) != 2 && chain[chain.size() - 2] == 0)
			fd = R.multiply(fd); // Apply a 180 deg x-rotation for atoms ''on the left'' of reference atom 1
	}

	fd.erase(fd.begin() + 3);

	return fd;
}

// Return second derivative of atom a1 (1-based) position
// with respect to Z-Matrix coordinates
// q2 = 1 (d), 2 (theta), 3 (phi) [first differentiation]
// q3 = 1 (d), 2 (theta), 3 (phi) [second differentiation]
// of, respectively, atoms a2 and a3 (1-based)
vectorOfDoubles molecule::getSecondDerivative(int a1, int q2, int a2, int q3, int a3)
{
	if ((a2 == 1) || (a2 == 2 && q2 != 1) || (a2 == 3 && q2 == 3))
	{
		std::cout << std::endl << "ERROR in molecule::getSecondDerivative(int a1, int q2, int a2, int q3, int a3) : no internal coordinate q2_a2 defined for q2 = " << q2 << " and a2 = " << a2 << std::endl << std::endl;
		exit(1);
	}
	else if ((a3 == 1) || (a3 == 2 && q3 != 1) || (a3 == 3 && q3 == 3))
	{
		std::cout << std::endl << "ERROR in molecule::getSecondDerivative(int a1, int q2, int a2, int q3, int a3) : no internal coordinate q3_a3 defined for q3 = " << q3 << " and a3 = " << a3 << std::endl << std::endl;
		exit(1);
	}

	int N1 = 0;
	for (int i = 0; i < nAtoms; ++i)
		N1 += atoms[i].getNFrag();

	int c1 = 4;
	int c2 = nAtoms - N1 + 2;

	matrix4D R(1.0,  0.0,  0.0,  0.0,
	    	   0.0, -1.0,  0.0,  0.0,
	   	   0.0,  0.0, -1.0,  0.0,
		   0.0,  0.0,  0.0,  1.0);

	int a1id = a1 - 1;
	int a2id = a2 - 1;
	int a3id = a3 - 1;

	vectorOfDoubles fd(4);
	vectorOfDoubles r0(4);

	matrix4D A(1.0, 0.0, 0.0, 0.0,
		   0.0, 1.0, 0.0, 0.0,
		   0.0, 0.0, 1.0, 0.0,
		   0.0, 0.0, 0.0, 1.0);

	fd[0] = 0.0; fd[1] = 0.0; fd[2] = 0.0; fd[3] = 0.0;
	r0[0] = 0.0; r0[1] = 0.0; r0[2] = 0.0; r0[3] = 1.0;

	switch (a1)
	{
		case 1:
		{
			if (a2 == 2 && q2 == 1 && a3 == 2 && q3 == 1)
			{
				A = atoms[atomsID[0]].getBmatSecondDerivative(1, 1);
				fd = A.multiply(r0);
			}
			break;
		}
		case 2:
		{
			// All derivatives are [0 0 0]
			// being atom 2 the reference atom
			break;
		}
		case 3:
		{
			if (a2 == 3 && a3 == 3)
			{
				A = atoms[atomsID[2]].getBmatSecondDerivative(q2, q3);
				fd = A.multiply(r0);
			}
			break;
		}
		case 4:
		{
			if (a2 == 3 && a3 == 3)
			{
				A = atoms[atomsID[3]].getBmat();
				A.multiply(atoms[atomsID[2]].getBmatSecondDerivative(q2, q3), A);
				fd = A.multiply(r0);
			}
			else if (a2 == 4 && a3 == 4)
			{
				A = atoms[atomsID[3]].getBmatSecondDerivative(q2, q3);
				A.multiply(atoms[atomsID[2]].getBmat(), A);
				fd = A.multiply(r0);
			}
			else if (a2 == 3 && a3 == 4)
			{
				A = atoms[atomsID[3]].getBmatFirstDerivative(q3);
				A.multiply(atoms[atomsID[2]].getBmatFirstDerivative(q2), A);
				fd = A.multiply(r0);
			}
			else if (a2 == 4 && a3 == 3)
			{
				A = atoms[atomsID[3]].getBmatFirstDerivative(q2);
				A.multiply(atoms[atomsID[2]].getBmatFirstDerivative(q3), A);
				fd = A.multiply(r0);
			}
			break;
		}
		default:
		{
			int a4id, a5id, count = 0;
			vectorOfIntegers chain = atoms[atomsID[a1id]].getChain();
			// Search if atom a2 is in the chain
			intit it = find(chain.begin(), chain.end(), a2id);
			if (it != chain.end())
			{
				if (a2 == a3) // If a3 = a2, then proceed
				{
					a4id = chain[0];
					while (a4id != a2id)
					{
						A.multiply(atoms[atomsID[a4id]].getBmat(), A);
						count++;
						a4id = chain[count];
					}
					A.multiply(atoms[atomsID[a2id]].getBmatSecondDerivative(q2, q3), A);
					if (a2id != 2)	// For a2id = 2 is intended atom 3. Thus, the premultiplication by A2 should not be applied (is a translation to atom 1)
					{
						a4id = atoms[atomsID[a2id]].getZMatAtom(1) - 1;
						A.multiply(atoms[atomsID[a4id]].getAmat(), A);
					}
					fd = A.multiply(r0);
				}
				else
				{
					count = 0;
					it = find(chain.begin(), chain.end(), a3id); // Search if a3 (!= a2) is in the chain
					if (it != chain.end()) // If both a2 and a3 are in the chain, proceed
					{
						for (int i = 0; i < chain.size(); ++i)
						{
							if (chain[i] == a2id)
							{
								A.multiply(atoms[atomsID[chain[i]]].getBmatFirstDerivative(q2), A);
								count++;
								if (count == 2)
								{
									if (a2id != 2)	// For a2id = 2 is intended atom 3. Thus, the premultiplication by A2 should not be applied (is a translation to atom 1)
									{
										a4id = atoms[atomsID[a2id]].getZMatAtom(1) - 1;
										A.multiply(atoms[atomsID[a4id]].getAmat(), A);
									}
									break;
								}
							}
							else if (chain[i] == a3id)
							{
								A.multiply(atoms[atomsID[chain[i]]].getBmatFirstDerivative(q3), A);
								count++;
								if (count == 2)
								{
									if (a3id != 2)	// For a3id = 2 is intended atom 3. Thus, the premultiplication by A2 should not be applied (is a translation to atom 1)
									{
										a4id = atoms[atomsID[a3id]].getZMatAtom(1) - 1;
										A.multiply(atoms[atomsID[a4id]].getAmat(), A);
									}
									break;
								}
							}
							else
								A.multiply(atoms[atomsID[chain[i]]].getBmat(), A);
						}
						fd = A.multiply(r0);
					}
				}
			}
		}
	}

	vectorOfIntegers chain;
	if (a1id >= c1 && a1id < c2)
	{
		chain = getAtom(a1id+1).getChain();
		if(atoms[atomsID[a1id]].getZMatAtom(1) != 2 && chain[chain.size() - 2] == 0)
			fd = R.multiply(fd); // Apply a 180 deg x-rotation for atoms ''on the left'' of reference atom 1
	}

	fd.erase(fd.begin() + 3);

	return fd;
}

// Build atoms chains from atom j to atoms 1, 2, 3, or 4
// and relative Aj matrices
void molecule::buildAmat(void)
{
	int a;
	vectorOfIntegers chain; // chain contains atoms ID in 0-based format
	matrix4D Amat;

	// Atom 1
	chain.clear();
	chain.push_back(0);
	chain.push_back(1);
	atoms[atomsID[0]].setChain(chain);
	atoms[atomsID[0]].setAmat(atoms[atomsID[0]].getBmat());

	// Atom 2
	chain.clear();
	chain.push_back(1);
	atoms[atomsID[1]].setChain(chain);
	atoms[atomsID[1]].setAmat(atoms[atomsID[1]].getBmat());

	// Atom 3
	chain.clear();
	chain.push_back(2);
	chain.push_back(1);
	atoms[atomsID[2]].setChain(chain);
	atoms[atomsID[2]].setAmat(atoms[atomsID[2]].getBmat());

	// Atom 4
	chain.clear();
	chain.push_back(3);
	chain.push_back(2);
	chain.push_back(1);
	atoms[atomsID[3]].setChain(chain);
	Amat = atoms[atomsID[3]].getBmat();
	Amat.multiply(atoms[atomsID[2]].getBmat(), Amat);
	atoms[atomsID[3]].setAmat(Amat);

	// Other atoms
	for (int i = 4; i < nAtoms; ++i)
	{
		chain.clear();
		chain.push_back(i);
		a = atoms[atomsID[i]].getZMatAtom(1);	
		while (a != 1 && a != 2 && a != 3 && a != 4)
		{
			a--;
			chain.push_back(a);
			a = atoms[atomsID[a]].getZMatAtom(1);
		}
		switch (a)
		{
			case 1:
			{
				chain.push_back(0);
				break;
			}
			case 2:
			{
				break;
			}
			case 3:
			{
				chain.push_back(2);
				break;
			}
			case 4:
			{
				chain.push_back(3);
				chain.push_back(2);
				break;
			}
		}
		chain.push_back(1);
		atoms[atomsID[i]].setChain(chain);

		Amat = atoms[atomsID[i]].getBmat();
		for (int j = 1; j < chain.size() - 1; ++j) // Avoid first (juist loaded in Amat) and last elements of chain
			Amat.multiply(atoms[atomsID[chain[j]]].getBmat(), Amat); 

		atoms[atomsID[i]].setAmat(Amat);
	}

	
	#ifdef DEBUG
	std::cout << "Atoms chains" << std::endl;
	for (int i = 0; i < nAtoms; ++i)
	{
		std::cout << i+1 << ") ";
		chain = atoms[atomsID[i]].getChain();
		for (int j = 0; j < chain.size(); ++j)
			std::cout << chain[j] + 1 << "\t";
		std::cout << std::endl;
	}
	std::cout << std::endl;
	#endif
	

	return;
}
	
// Find the min ID atom bonded to given atom
// NB: input atom ID is 0-based, while output atom ID is 1-based
// a: 0-based id
// b: 1-based id's
int molecule::minBond(int a, vectorOfIntegers b)
{
	int check;
	int nb = b.size();
        // SR 15-02-21 BEGIN //
        // Possible bug fix. m was initialized here as atoms.size() + 1. If the following conditions were not met, the procedure returned atomsID[m-1]+1, but atomsID is 0-based and atomsID[atoms.size()] is out of range.
	// int m = atoms.size() + 1;
	int m = atoms.size();
        // SR 15-02-21 END   //
	for (int i = 0; i < atoms[a].getNBonds(); ++i)
	{
		check = 1;
		for (int j = 0; j < nb; ++j)
		{
			if (atoms[a].getBond(i) == b[j])
				check = 0;
		}
		//if (atoms[a].id == 2 && nb == 2 && atoms[atoms[a].getBond(i)-1].id == 1) check = 0; // This is to avoid the X-2-1-... Z-matrix enty for atoms X directly connected to 2. Instead, the X-2-3-4 Z-Matrix is created
		m = check ? MIN(m, atoms[atoms[a].getBond(i)-1].id) : m;
	}
	return (atomsID[m-1]+1);
}

// This method permutes atoms so that they
// are ordered following the convention for
// the Z-Matrix construction
//
// NB: ID's in atomsID are 0-based
void molecule::permute()
{
	int N1 = 0;
	for (int i = 0; i < nAtoms; ++i)
		N1 += atoms[i].getNFrag();

	int c1 = 4;
	int c2 = nAtoms - N1 + 2;

	for (int i = 0; i < nAtoms; ++i)
	{
		atomsID[i] = -1;
		atoms[i].id = -1;
		atoms[i].activeID = -1;
	}

	for (int i = 0; i < nActiveAtoms; ++i)
		atomsActiveID[i] = -1;

	atoms[dA1 - 1].id = 1; atoms[dA1 - 1].activeID = 1;
	atoms[dA2 - 1].id = 2; atoms[dA2 - 1].activeID = 2;
	atoms[dA3 - 1].id = 3; atoms[dA3 - 1].activeID = 3;
	atoms[dA4 - 1].id = 4; atoms[dA4 - 1].activeID = 4;

	int stop = 0;
	int aid, c0 = 5, activeC0 = 5;
	vectorOfIntegers bonds;

	// Fragment 0
	while (!stop)
	{
		stop = 1;

		// search for next atom in fragment 0
		for (int i = 0; i < nAtoms; ++i)
		{
			if (atoms[i].id > -1 && atoms[i].getNFrag() == 0 && atomsID[i] == -1)
			{
				aid = i;
				atomsID[i] = 0;
				stop = 0;
				break;
			}
		}

		// Assign ID to atoms bonded to atom aid
		bonds = atoms[aid].getBonds();
		for (int i = 0; i < bonds.size(); ++i)
		{
			if (atoms[bonds[i]-1].id == -1)
			{
				atoms[bonds[i]-1].id = c0;
				c0++;
				if (atoms[bonds[i]-1].isActive())
				{
					atoms[bonds[i]-1].activeID = activeC0;
					activeC0++;
				}
			}
		}
	}

	// Fragment 1
	stop = 0;
	while (!stop)
	{
		stop = 1;

		// search for next atom in fragment 0
		for (int i = 0; i < nAtoms; ++i)
		{
			if (atoms[i].id > -1 && atoms[i].getNFrag() == 1 && atomsID[i] == -1)
			{
				aid = i;
				atomsID[i] = 0;
				stop = 0;
				break;
			}
		}

		// Assign ID to atoms bonded to atom aid
		bonds = atoms[aid].getBonds();
		for (int i = 0; i < bonds.size(); ++i)
		{
			if (atoms[bonds[i]-1].id == -1)
			{
				atoms[bonds[i]-1].id = c0;
				c0++;
				if (atoms[bonds[i]-1].isActive())
				{
					atoms[bonds[i]-1].activeID = activeC0;
					activeC0++;
				}
			}
		}
	}

	// Fill atomsID array

	int ai = 0;
	for (int i = 0; i < nAtoms; ++i)
	{
		atomsID[atoms[i].id-1] = i;
		if (atoms[i].isActive())
		{
			atomsActiveID[atoms[i].activeID-1] = ai;
			ai++;
		}
	}
	
	return;
}

// Calculate internal coordinates among 2, 3 or 4 atoms
// NB: atoms ID are 1-based
double molecule::getBondLength(int a1, int a2)
{
	double dx = atoms[a2-1].x() - atoms[a1-1].x();
	double dy = atoms[a2-1].y() - atoms[a1-1].y();
	double dz = atoms[a2-1].z() - atoms[a1-1].z();

	return sqrt(dx * dx + dy * dy + dz * dz);
}

double molecule::getBondAngle(int a1, int a2, int a3)
{
	double v1x = atoms[a2-1].x() - atoms[a1-1].x();
	double v1y = atoms[a2-1].y() - atoms[a1-1].y();
	double v1z = atoms[a2-1].z() - atoms[a1-1].z();

	double v2x = atoms[a2-1].x() - atoms[a3-1].x();
	double v2y = atoms[a2-1].y() - atoms[a3-1].y();
	double v2z = atoms[a2-1].z() - atoms[a3-1].z();

	double n1 = sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
	double n2 = sqrt(v2x * v2x + v2y * v2y + v2z * v2z);

	double t = acos((v1x * v2x + v1y * v2y + v1z * v2z) / (n1 * n2));

	// Wrap angle in (-180, 180] range
	
	t = t - 2.0 * M_PI * floor((t + M_PI) / (2.0 * M_PI));

	return t;
}

double molecule::getDihedralAngle(int a1, int a2, int a3, int a4)
{
	double b1x = atoms[a2-1].x() - atoms[a1-1].x();
	double b1y = atoms[a2-1].y() - atoms[a1-1].y();
	double b1z = atoms[a2-1].z() - atoms[a1-1].z();

	double b2x = atoms[a3-1].x() - atoms[a2-1].x();
	double b2y = atoms[a3-1].y() - atoms[a2-1].y();
	double b2z = atoms[a3-1].z() - atoms[a2-1].z();

	double b3x = atoms[a4-1].x() - atoms[a3-1].x();
	double b3y = atoms[a4-1].y() - atoms[a3-1].y();
	double b3z = atoms[a4-1].z() - atoms[a3-1].z();

	double v1x = b1y * b2z - b1z * b2y;
	double v1y = b1z * b2x - b1x * b2z;
	double v1z = b1x * b2y - b1y * b2x;

	double v2x = b2y * b3z - b2z * b3y;
	double v2y = b2z * b3x - b2x * b3z;
	double v2z = b2x * b3y - b2y * b3x;

        double ux = v1y * v2z - v1z * v2y;
        double uy = v1z * v2x - v1x * v2z;
        double uz = v1x * v2y - v1y * v2x;

        double nb2 = sqrt(b2x * b2x + b2y * b2y + b2z * b2z);

        double w = v1x * v2x + v1y * v2y + v1z * v2z;
        double u = (ux * b2x + uy * b2y + uz * b2z) / nb2;
        double t = atan2(u, w);

        // Wrap angle in (-180, 180] range

	t = t - 2.0 * M_PI * floor((t + M_PI) / (2.0 * M_PI));

	return t;
}

// Fragment the molecule in two pieces

void molecule::fragmentMolecule(void)
{
	int i, j, k, l, torsCounter;
	int at1, at2, target, *scatter;
	int *atmxfrg;
	int n_dihedrals_old, nfragments_old, fr1, fr2;
	char string[255], stmp[10], end;
	double dtmp;
	FILE *f1;

	/******************/
	/* Make fragments */
	/******************/

	fragments = (frag *)calloc(nfragments, sizeof(frag));
	fragments[0].isref = 1; fragments[0].next = -1; fragments[0].nt = -1;
	fragments[1].isref = 0; fragments[1].next = 1;  fragments[1].nt =  0;
	
	at1 = dA2; at2 = dA3;
	at1--; at2--;

	k = 0;
	j = atoms[at1].getNBonds();
	for (i = 0; i < j; i++)
	{
		if (atoms[at1].getBond(i) == at2 + 1)
		{
			#ifdef _DEBUGGING_
			printf("\nOk, atoms %d and %d are connected so we can proceed!",at1+1,at2+1);
			#endif
			k = 1;
			break;
		}
	}
	if (!k)
	{
		printf("\n\nERROR : atoms %d and %d are not bonded!\n\n",++at1,++at2);
		exit(0);
	}

	scatter = (int *)calloc(nAtoms, sizeof(int));
	for (i = 0; i < nAtoms; i++) scatter[i] = -1;

	atoms[at1].setNFrag(1);
	atoms[at2].setNFrag(2);
	scatter[at1] = scatter[at2] = 0;

	for (i = 0; i < atoms[at1].getNBonds(); i++)
	{
		j = atoms[at1].getBond(i) - 1;
		if (j != at2) scatter[j] = 1;
	}

	for (i = 0; i <atoms[at2].getNBonds(); i++)
	{
		j = atoms[at2].getBond(i) - 1;
		if (j != at1) scatter[j] = 2;
	}

	do {
		k = 0;
		for (i = 0; i < nAtoms; i++)
		{
			if (scatter[i] > 0)
			{
				k=1;
				atoms[i].setNFrag(scatter[i]);
				for (j = 0; j < atoms[i].getNBonds(); j++)
				{
					if (scatter[atoms[i].getBond(j) - 1] == -1)
					{
						scatter[atoms[i].getBond(j) - 1] = scatter[i];
					}
				}
				scatter[i] = 0;
			}
		}
  	} while(k);

 	j = 0;
 	for (i = 0; i < nAtoms; i++)
	{
 		if (atoms[i].getNFrag() == 1) j++;
	}


// CODE TO BE ADAPTED IN CASE OF MORE DIHEDRAL ANGLES
//	n_dihedrals_old = 1;
//	nfragments_old = nfragments;
//	l = 0;
//	torsCounter = 1;
//	do {
//	
//	    at1 = pep_bonds[torsCounter*2+0];
//	    at2 = pep_bonds[torsCounter*2+1];
//	    if(!at1 || !at2) break;
//	    at1--; at2--;
//	    j=atoms[at1].nbonds;
//	    k=0;
//	    for (i=0;i<j;i++){
//	      if (*(atoms[at1].bonds+i)==at2+1){
//	#ifdef _DEBUGGING_
//		printf("\nOk, atoms %d and %d are connected so we can proceed!",at1+1,at2+1);
//	#endif
//		k=1;
//		break;
//	      }
//	    }
//	    if (!k){
//	      printf("\n\nERROR : atoms %d and %d are not bonded!\n\n",++at1,++at2);
//	      exit(0);
//	    }
//	
//	    nfragments++;
//	    fragments=(frag *)realloc(fragments,(nfragments+1)*sizeof(frag));
//	    fr2=atoms[at1].nfrag;
//	
//	    j=atoms[at1].nfrag-1;
//	    for (i=0;i<natoms;i++){
//	      if (atoms[i].nfrag-1!=j) scatter[i]=-2;
//	      else if (i==at1) scatter[i]=-3;
//	      else if (i==at2) scatter[i]=-4;
//	      else scatter[i]=-1;
//	    }
//	  
//	    for (i=0;i<atoms[at1].nbonds;i++){
//	      j=*(atoms[at1].bonds+i)-1;
//	      if (j!=at2) scatter[j]=-3;
//	      if(atoms[at1].nfrag != atoms[j].nfrag) scatter[j]=-2;
//	    }
//	    for (i=0;i<atoms[at2].nbonds;i++){
//	      j=*(atoms[at2].bonds+i)-1;
//	      if (j!=at1) scatter[j]=-4;
//	      if (atoms[at2].nfrag != atoms[j].nfrag) scatter[j]=-2;
//	    }
//	
//	    fr2=atoms[at1].nfrag;
//	    if (!fragments[fr2-1].isref){
//	      i=fragments[fr2-1].nt;
//	      target=pep_bonds[i*2+1]-1;
//	    }
//	    else{
//	      for (i=0;i<torsCounter;i++){
//		if (atoms[pep_bonds[i*2+0]-1].nfrag==fr2 & pep_bonds[i*2+0]!=at1+1){
//		  target=pep_bonds[i*2+0]-1;
//		  break;
//		}
//	      }
//	    }
//	
//	    do{
//	      k=0;
//	      for (i=0;i<natoms;i++){
//		if (scatter[i]<-2){
//		  k=1;
//		  atoms[i].nfrag=scatter[i];
//		  for (j=0;j<atoms[i].nbonds;j++){
//		    if (scatter[*(atoms[i].bonds+j)-1]==-1){
//		      scatter[*(atoms[i].bonds+j)-1]=scatter[i];
//		    }
//	
//		    if (*(atoms[i].bonds+j)-1==target)
//		      fr1=scatter[i];
//		    
//		  }
//		  scatter[i]=0;
//		}
//	      }
//	    }while(k);
//	  
//	    k=(fr1==-3 ? -4 : -3);
//	    for (i=0;i<natoms;i++){
//	      if (atoms[i].nfrag==fr1) atoms[i].nfrag=fr2;
//	      else if (atoms[i].nfrag==k) atoms[i].nfrag=nfragments;
//	    }
//	  
//	    fragments[nfragments-1].next=fr2;
//	    fragments[nfragments-1].nt=torsCounter;///////-1;
//	
//	    if (atoms[at1].nfrag==nfragments){
//	      i=pep_bonds[(torsCounter-1)*2+0];
//	      pep_bonds[(torsCounter-1)*2+0]=pep_bonds[(torsCounter-1)*2+1];
//	      pep_bonds[(torsCounter-1)*2+1]=i;
//	    }
//	
//	    n_dihedrals_old=torsCounter;
//	    nfragments_old=nfragments;
//		
//		torsCounter++;
//	  
//	  }while(!l && torsCounter<n_dihedrals);

 _JUMP:

	for (i=0; i < nAtoms; i++)
		atoms[i].setNFrag(atoms[i].getNFrag() - 1);
	for (i = 0; i < nfragments; i++)
		fragments[i].next--;
  
 	/* //DB//
	printf("\nDISTRIBUTION OF ATOMS OVER FRAGMENTS\n\n");
	printf("Atom\tFragment\n\n");
	for (i = 0; i < nAtoms; i++)
	{
		printf("%4d\t%8d\n", i + 1, atoms[i].getNFrag());
  	}
 	printf("\n\nMOLECULAR TOPOLOGY\n\n");
 	printf("Fragment\tIsRef\tNext Frag.\tTorsional\n\n");
 	for (i = 0; i < nfragments; i++)
	{
 		printf("%8d\t%5d\t%10d\t%9d\n",i, fragments[i].isref, fragments[i].next, fragments[i].nt);
	}
	printf("\n\nROTATEBLE BONDS ARE DEFINED BY\n\n");
	printf("Torsional\tAtom1 ID\tAtom2 ID\n\n");
	for (i = 0; i < n_dihedrals; i++)
		printf("%9d\t%8d\t%8d\n", i, dA2, dA3);
	std::cout << std::endl;
	
	 OUTPUT TO FILE
	#ifdef _LINUX_
	sprintf(string, "%s/%s_fragments_%d", path, project, (imol+1));
	#else
	sprintf(string,"%s\\%s_fragments_%d",path,project,(imol+1));
	#endif
	f1=fopen(string,"w");
	for ( i = 0; i < nAtoms; i++)
	{
		fprintf(f1, "%d\n", atoms[i].getNFrag());
	}
	fclose(f1);
	*/

  return;

}


////////////////////////////////////////////////////////////////

// TO CANCEL ?

//	void single_fragment(int imol)
//	{
//
//		int i_index;
//		char str[10000];
//		FILE *f1;
//
//	#ifdef _LINUX_
//		sprintf(str,"%s/%s_fragments_%d",path,project,(imol+1));
//	#else
//		sprintf(str,"%s\\%s_fragments_%d",path,project,(imol+1));
//	#endif
//		nfragments=1;
//		fragments=(frag *)calloc(1,sizeof(frag));
//		fragments[0].isref=1;
//		fragments[0].next=fragments[0].nt=-1;
//		fragments[0].natoms = natoms;
//		fragments[0].atoms = (int*)calloc(natoms,sizeof(int));
//		f1=fopen(str,"w");
//		for (i_index=0;i_index<natoms;i_index++){
//			atoms[i_index].nfrag=0;
//			fragments[0].atoms[i_index] = i_index;
//			fprintf(f1,"%d\n",atoms[i_index].nfrag);
//		}
//		fclose(f1);
//		
//		// Set the fragment as terminal
//		nterm[0] = nterm[1] = 0;
//		// Set the ordering array
//		ordering = (int*)calloc(1,sizeof(int));
//		ordering[0] = 0;
//		return;
//	}

/***************************/
/* Scale atoms coordinates */
/***************************/

void molecule::scaleAtomsPositions(double scale)
{
	for (int i = 0; i < nAtoms; ++i)
	{
		atoms[i].x(scale * atoms[i].x());
		atoms[i].y(scale * atoms[i].y());
		atoms[i].z(scale * atoms[i].z());
	}
	return;
}

/****************************************/
/* Check if refence atoms have been set */
/****************************************/

int molecule::areRefAtomsSet(void)
{
	return referenceAtomsAreSet;
}

int molecule::hasZMatrix(void)
{
	return isZMatrixBuilt;
}

/*********************************************************************************/
/* Build and store Cartesian coordinates, and first and second derivatives in MF */
/*********************************************************************************/
inline int molecule::qidx(int qid)
{
	switch(qid)
	{
		case 3:
			return 0;
		case 6:
			return 1;
		case 7:
			return 2;
		default:
			return (qid - 6);
	}
}

void molecule::calculateDerivativesOfCM(void)
{
	int nx = 3 * nAtoms;
	int nq = 3 * nActiveAtoms - 6;
	atom at;

	// 1. Store the position of the CM and
	//    the matrix that provides the rotation
	//    BF --> CF (i.e., MF)
	//    where BF is obtained from AF by means
	//    of the translation of the CM.
	calculateCenterOfMass(); // this is the position of the CM

	calculateInertiaTensorInCM();
	matrix3D inert = getInertiaTensor();
	
 	diag.setMaxStep(500);
	diag.setMatrix(inert);
	diag.diagonalize();
	diag.reorder();
	matrix3D Emat = diag.getEigenVectors3D();
	Emat.transpose();
	rotationToMF = Emat; // this is the rotation matrix


	// 2. First derivatives of CM

	double M = 0.0, atMass[nAtoms];
	dRCM = vectorOfDoubles(nq * 3, 0.0);
	vectorOfDoubles d(3);

	// q = d2
	dRCM[0 * 3 + 0] = 0.0;
	dRCM[0 * 3 + 1] = 0.0;
	dRCM[0 * 3 + 2] = 0.0;
	for (int i = 1; i <= nAtoms; ++i)
	{ 
		at = getAtom(i);
		atMass[i] = at.getMass();
		M += atMass[i];
		d = getFirstDerivative(i, 1, 2);
		dRCM[0 * 3 + 0] += atMass[i] * d[0];
		dRCM[0 * 3 + 1] += atMass[i] * d[1];
		dRCM[0 * 3 + 2] += atMass[i] * d[2];
	}
	M = 1.0 / M;
	dRCM[0 * 3 + 0] *= M;
	dRCM[0 * 3 + 1] *= M;
	dRCM[0 * 3 + 2] *= M;

	// q = d3
	dRCM[1 * 3 + 0] = 0.0;
	dRCM[1 * 3 + 1] = 0.0;
	dRCM[1 * 3 + 2] = 0.0;
	for (int i = 1; i <= nAtoms; ++i)
	{ 
		d = getFirstDerivative(i, 1, 3);
		dRCM[1 * 3 + 0] += atMass[i] * d[0];
		dRCM[1 * 3 + 1] += atMass[i] * d[1];
		dRCM[1 * 3 + 2] += atMass[i] * d[2];
	}
	dRCM[1 * 3 + 0] *= M;
	dRCM[1 * 3 + 1] *= M;
	dRCM[1 * 3 + 2] *= M;

	// q = theta3
	dRCM[2 * 3 + 0] = 0.0;
	dRCM[2 * 3 + 1] = 0.0;
	dRCM[2 * 3 + 2] = 0.0;
	for (int i = 1; i <= nAtoms; ++i)
	{ 
		d = getFirstDerivative(i, 2, 3);
		dRCM[2 * 3 + 0] += atMass[i] * d[0];
		dRCM[2 * 3 + 1] += atMass[i] * d[1];
		dRCM[2 * 3 + 2] += atMass[i] * d[2];
	}
	dRCM[2 * 3 + 0] *= M;
	dRCM[2 * 3 + 1] *= M;
	dRCM[2 * 3 + 2] *= M;

	// All other internal coordinates

	int idj;
	atom *activeAtoms = new atom[nActiveAtoms];
	getActiveAtoms(activeAtoms);

	for (int j = 3; j < nq; j += 3)
	{
		dRCM[j * 3 + 0] = 0.0;
		dRCM[j * 3 + 1] = 0.0;
		dRCM[j * 3 + 2] = 0.0;

		dRCM[(j + 1) * 3 + 0] = 0.0;
		dRCM[(j + 1) * 3 + 1] = 0.0;
		dRCM[(j + 1) * 3 + 2] = 0.0;

		dRCM[(j + 2) * 3 + 0] = 0.0;
		dRCM[(j + 2) * 3 + 1] = 0.0;
		dRCM[(j + 2) * 3 + 2] = 0.0;

		idj = activeAtoms[3 + j / 3 - 1].getID();

		for (int i = 1; i <= nAtoms; ++i)
		{ 
			d = getFirstDerivative(i, 1, idj);
			dRCM[j * 3 + 0] += atMass[i] * d[0];
			dRCM[j * 3 + 1] += atMass[i] * d[1];
			dRCM[j * 3 + 2] += atMass[i] * d[2];

			d = getFirstDerivative(i, 2, idj);
			dRCM[(j + 1) * 3 + 0] += atMass[i] * d[0];
			dRCM[(j + 1) * 3 + 1] += atMass[i] * d[1];
			dRCM[(j + 1) * 3 + 2] += atMass[i] * d[2];

			d = getFirstDerivative(i, 3, idj);
			dRCM[(j + 2) * 3 + 0] += atMass[i] * d[0];
			dRCM[(j + 2) * 3 + 1] += atMass[i] * d[1];
			dRCM[(j + 2) * 3 + 2] += atMass[i] * d[2];
		}

		dRCM[j * 3 + 0] *= M;
		dRCM[j * 3 + 1] *= M;
		dRCM[j * 3 + 2] *= M;

		dRCM[(j + 1) * 3 + 0] *= M;
		dRCM[(j + 1) * 3 + 1] *= M;
		dRCM[(j + 1) * 3 + 2] *= M;

		dRCM[(j + 2) * 3 + 0] *= M;
		dRCM[(j + 2) * 3 + 1] *= M;
		dRCM[(j + 2) * 3 + 2] *= M;
	}

	// 3. Second derivatives of CM

	d2RCM = vectorOfDoubles(nq * nq * 3, 0.0);

	/* 
	 * The following array contains the id's of q coordinates
	 *
	 * Qid = atID_0 * 3 + qType
	 *
	 * atID_0 : 0-based atom ID (i.e., line in the Z-Matrix)
	 * qType = 0 (d), 1 (theta), 2 (phi)
	 */
	int *Qid = new int[nq];
	Qid[0] = 3; Qid[1] = 6; Qid[2] = 7;
	for (int i = 3; i < nq; i += 3)
	{
		Qid[i + 0] = i * 3 + 0; 
		Qid[i + 1] = i * 3 + 1; 
		Qid[i + 2] = i * 3 + 2; 
	}

	int q1, a1, q2, a2;
	vectorOfIntegers chain;

	for (int i = 1; i <= nAtoms; ++i)
	{
		at = getAtom(i);
		atMass[i] = at.getMass() * M;  // is the ratio: mi / M
		chain = at.getChain();
		for (int j1 = 0; j1 < chain.size(); ++j1)
		{
			if (getAtom(chain[j1] + 1).isActive())
			{
				for (int j2 = 0; j2 < chain.size(); ++j2)
				{
					if (getAtom(chain[j2] + 1).isActive())
					{
						q1 = qidx(getAtom(chain[j1] + 1).activeID * 3 + 0);
						q2 = qidx(getAtom(chain[j2] + 1).activeID * 3 + 0);

						// dj1 - dj2
						if (chain[j1] > 1 && chain[j2] > 1)
						{
							d = getSecondDerivative(i, 1, chain[j1], 1, chain[j2]);
							d2RCM[(q1 * nq + q2) * 3 + 0] += atMass[i] * d[0];
							d2RCM[(q1 * nq + q2) * 3 + 1] += atMass[i] * d[1];
							d2RCM[(q1 * nq + q2) * 3 + 2] += atMass[i] * d[2];
						}
						// dj1 - thetaj2
						if (chain[j1] > 1 && chain[j2] > 2)
						{
							d = getSecondDerivative(i, 1, chain[j1], 2, chain[j2]);
							d2RCM[(q1 * nq + q2) * 3 + 0] += atMass[i] * d[0];
							d2RCM[(q1 * nq + q2) * 3 + 1] += atMass[i] * d[1];
							d2RCM[(q1 * nq + q2) * 3 + 2] += atMass[i] * d[2];
						}
						// dj1 - phij2
						if (chain[j1] > 1 && chain[j2] > 3)
						{
							d = getSecondDerivative(i, 1, chain[j1], 3, chain[j2]);
							d2RCM[(q1 * nq + q2) * 3 + 0] += atMass[i] * d[0];
							d2RCM[(q1 * nq + q2) * 3 + 1] += atMass[i] * d[1];
							d2RCM[(q1 * nq + q2) * 3 + 2] += atMass[i] * d[2];
						}
						// thetaj1 - dj2
						if (chain[j1] > 2 && chain[j2] > 1)
						{
							d = getSecondDerivative(i, 2, chain[j1], 1, chain[j2]);
							d2RCM[(q1 * nq + q2) * 3 + 0] += atMass[i] * d[0];
							d2RCM[(q1 * nq + q2) * 3 + 1] += atMass[i] * d[1];
							d2RCM[(q1 * nq + q2) * 3 + 2] += atMass[i] * d[2];
						}
						// thetaj1 - thetaj2
						if (chain[j1] > 2 && chain[j2] > 2)
						{
							d = getSecondDerivative(i, 2, chain[j1], 2, chain[j2]);
							d2RCM[(q1 * nq + q2) * 3 + 0] += atMass[i] * d[0];
							d2RCM[(q1 * nq + q2) * 3 + 1] += atMass[i] * d[1];
							d2RCM[(q1 * nq + q2) * 3 + 2] += atMass[i] * d[2];
						}
						// thetaj1 - phij2
						if (chain[j1] > 2 && chain[j2] > 3)
						{
							d = getSecondDerivative(i, 2, chain[j1], 3, chain[j2]);
							d2RCM[(q1 * nq + q2) * 3 + 0] += atMass[i] * d[0];
							d2RCM[(q1 * nq + q2) * 3 + 1] += atMass[i] * d[1];
							d2RCM[(q1 * nq + q2) * 3 + 2] += atMass[i] * d[2];
						}
						// phi1 - dj2
						if (chain[j1] > 3 && chain[j2] > 1)
						{
							d = getSecondDerivative(i, 3, chain[j1], 1, chain[j2]);
							d2RCM[(q1 * nq + q2) * 3 + 0] += atMass[i] * d[0];
							d2RCM[(q1 * nq + q2) * 3 + 1] += atMass[i] * d[1];
							d2RCM[(q1 * nq + q2) * 3 + 2] += atMass[i] * d[2];
						}
						// phi1 - thetaj2
						if (chain[j1] > 3 && chain[j2] > 2)
						{
							d = getSecondDerivative(i, 3, chain[j1], 2, chain[j2]);
							d2RCM[(q1 * nq + q2) * 3 + 0] += atMass[i] * d[0];
							d2RCM[(q1 * nq + q2) * 3 + 1] += atMass[i] * d[1];
							d2RCM[(q1 * nq + q2) * 3 + 2] += atMass[i] * d[2];
						}
						// phij1 - phij2
						if (chain[j1] > 3 && chain[j2] > 3)
						{
							d = getSecondDerivative(i, 3, chain[j1], 3, chain[j2]);
							d2RCM[(q1 * nq + q2) * 3 + 0] += atMass[i] * d[0];
							d2RCM[(q1 * nq + q2) * 3 + 1] += atMass[i] * d[1];
							d2RCM[(q1 * nq + q2) * 3 + 2] += atMass[i] * d[2];
						}
					}
				}
			}
		}
	}
	return;
}

// RETURN THE CARTESIAN COORDINATE OF ATOM atID EXPRESSED IN MF
// atID : 1-based atom ID
vectorOfDoubles molecule::getAtomPositionInMF(int atID)
{
	vectorOfDoubles r0(3), r(3);
	r0 = getAtomPosition(atID);
	r[0] = rotationToMF.xx * (r0[0] - centerOfMass[0]) + rotationToMF.xy * (r0[1] - centerOfMass[1]) + rotationToMF.xz * (r0[2] - centerOfMass[2]);
	r[1] = rotationToMF.yx * (r0[0] - centerOfMass[0]) + rotationToMF.yy * (r0[1] - centerOfMass[1]) + rotationToMF.yz * (r0[2] - centerOfMass[2]);
	r[2] = rotationToMF.zx * (r0[0] - centerOfMass[0]) + rotationToMF.zy * (r0[1] - centerOfMass[1]) + rotationToMF.zz * (r0[2] - centerOfMass[2]);
	return r;
}

// RETURN THE FIRST DERIVATIVE, EXPRESSED IN MF, OF ATOM a1 (1-based)
// WITH RESPECT TO Z-Matrix COORDINATE
// q = 1 (d), 2 (theta), 3 (phi)
// OF ATOM a2 (1-based)
vectorOfDoubles molecule::getFirstDerivativeInMF(int a1, int q, int a2)
{
	// THESE FORMULAS ARE BASED ON THE FACT THAT THE AF --> MF ROTO-TRANSLATION IS CONSTANT AT FIXED MOLECULAR CONFIGURATION
	// THUS, ONLY ROTATION IS APPLIED
	vectorOfDoubles d0(3), d(3);
	d0 = getFirstDerivative(a1, q, a2);
	d[0] = rotationToMF.xx * d0[0] + rotationToMF.xy * d0[1] + rotationToMF.xz * d0[2];
	d[1] = rotationToMF.yx * d0[0] + rotationToMF.yy * d0[1] + rotationToMF.yz * d0[2];
	d[2] = rotationToMF.zx * d0[0] + rotationToMF.zy * d0[1] + rotationToMF.zz * d0[2];
	return d;
}

// Return second derivative, EXPRESSED IN MF, of atom a1 (1-based) position
// with respect to Z-Matrix coordinates
// q2 = 1 (d), 2 (theta), 3 (phi) [first differentiation]
// q3 = 1 (d), 2 (theta), 3 (phi) [second differentiation]
// of, respectively, atoms a2 and a3 (1-based)
vectorOfDoubles molecule::getSecondDerivativeInMF(int a1, int q2, int a2, int q3, int a3)
{
	// THESE FORMULAS ARE BASED ON THE FACT THAT THE AF --> MF ROTO-TRANSLATION IS CONSTANT AT FIXED MOLECULAR CONFIGURATION
	// THUS, ONLY ROTATION IS APPLIED
	vectorOfDoubles d20(3), d2(3);
	d20 = getSecondDerivative(a1, q2, a2, q3, a3);
	d2[0] = rotationToMF.xx * d20[0] + rotationToMF.xy * d20[1] + rotationToMF.xz * d20[2];
	d2[1] = rotationToMF.yx * d20[0] + rotationToMF.yy * d20[1] + rotationToMF.yz * d20[2];
	d2[2] = rotationToMF.zx * d20[0] + rotationToMF.zy * d20[1] + rotationToMF.zz * d20[2];
	return d2;
}

// THE FOLLOWING METHODS ALLOW ONE TO OBTAIN THE COORDINATES AND DERIVATIVES
// OF ATOMS IN A FRAME WITH A CUSTOM ORIENTATION WITH RESPECT TO MF. THIS
// CAN BE OBTAINED BY DEFINING ANY SET OF EULER ANGLES OR USING 3 ATOMS TO
// BUILD A NEW REFERENCE SYSTEM
void molecule::setCustomOrientation(double al, double be, double ga)
{
	double ca = cos(al), sa = sin(al);
	double cb = cos(be), sb = sin(be);
	double cg = cos(ga), sg = sin(ga);

	customRot.xx = ca * cb * cg - sa * sg;
	customRot.xy = sa * cb * cg + ca * sg;
	customRot.xz = -sb * cg;
	customRot.yx = -ca * cb * sg - sa * cg;
	customRot.yy = -sa * cb * sg + ca * cg;
	customRot.yz = sb * sg;
	customRot.zx = ca * sb;
	customRot.zy = sa * sb;
	customRot.zz = cb;
	return;
}
void molecule::setCustomOrientation(int a1, int a2, int a3)
{
	vectorOfDoubles rA1 = getAtomPosition(a1);
	vectorOfDoubles rA2 = getAtomPosition(a2);
	vectorOfDoubles rA3 = getAtomPosition(a3);

	double n;
	vectorOfDoubles z(3), v(3), x(3), y(3);

	z[0] = rA1[0] - rA2[0];
	z[1] = rA1[1] - rA2[1];
	z[2] = rA1[2] - rA2[2];
	n = 1.0 / sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
	z[0] *= n;
	z[1] *= n;
	z[2] *= n;

	v[0] = rA3[0] - rA2[0];
	v[1] = rA3[1] - rA2[1];
	v[2] = rA3[2] - rA2[2];

	y[0] = v[1] * z[2] - v[2] * z[1];
	y[1] = v[2] * z[0] - v[0] * z[2];
	y[2] = v[0] * z[1] - v[1] * z[0];
	n = 1.0 / sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
	y[0] *= n;
	y[1] *= n;
	y[2] *= n;

	x[0] = y[1] * z[2] - y[2] * z[1];
	x[1] = y[2] * z[0] - y[0] * z[2];
	x[2] = y[0] * z[1] - y[1] * z[0];

	customRot.xx = x[0];
	customRot.xy = x[1];
	customRot.xz = x[2];
	customRot.yx = y[0];
	customRot.yy = y[1];
	customRot.yz = y[2];
	customRot.zx = z[0];
	customRot.zy = z[1];
	customRot.zz = z[2];

	return;
}

void molecule::setCustomOrientation(matrix3D m)
{
	customRot.xx = m.xx;
	customRot.xy = m.xy;
	customRot.xz = m.xz;
	customRot.yx = m.yx;
	customRot.yy = m.yy;
	customRot.yz = m.yz;
	customRot.zx = m.zx;
	customRot.zy = m.zy;
	customRot.zz = m.zz;
	return;
}

void molecule::setCustomCenter(vectorOfDoubles T)
{
	customCenter[0] = T[0];
	customCenter[1] = T[1];
	customCenter[2] = T[2];
	return;
}

void molecule::setCustomCenter(int atomID)
{
	customCenter = getAtomPosition(atomID);
	return;
}

vectorOfDoubles molecule::getAtomPositionInCustomFrame(int atID)
{
	vectorOfDoubles r0(3), r(3);
	r0 = getAtomPosition(atID);
	r[0] = customRot.xx * (r0[0] - customCenter[0]) + customRot.xy * (r0[1] - customCenter[1]) + customRot.xz * (r0[2] - customCenter[2]); 
	r[1] = customRot.yx * (r0[0] - customCenter[0]) + customRot.yy * (r0[1] - customCenter[1]) + customRot.yz * (r0[2] - customCenter[2]);
	r[2] = customRot.zx * (r0[0] - customCenter[0]) + customRot.zy * (r0[1] - customCenter[1]) + customRot.zz * (r0[2] - customCenter[2]);
	return r;
}

vectorOfDoubles molecule::getFirstDerivativeInCustomFrame(int a1, int q, int a2)
{
	vectorOfDoubles d0(3), d(3);
	d0 = getFirstDerivative(a1, q, a2);
	d[0] = customRot.xx * d0[0] + customRot.xy * d0[1] + customRot.xz * d0[2];
	d[1] = customRot.yx * d0[0] + customRot.yy * d0[1] + customRot.yz * d0[2];
	d[2] = customRot.zx * d0[0] + customRot.zy * d0[1] + customRot.zz * d0[2];
	return d;
}

vectorOfDoubles molecule::getSecondDerivativeInCustomFrame(int a1, int q2, int a2, int q3, int a3)
{
	vectorOfDoubles d20(3), d2(3);
	d20 = getSecondDerivative(a1, q2, a2, q3, a3);
	d2[0] = customRot.xx * d20[0] + customRot.xy * d20[1] + customRot.xz * d20[2];
	d2[1] = customRot.yx * d20[0] + customRot.yy * d20[1] + customRot.yz * d20[2];
	d2[2] = customRot.zx * d20[0] + customRot.zy * d20[1] + customRot.zz * d20[2];
	return d2;
}

