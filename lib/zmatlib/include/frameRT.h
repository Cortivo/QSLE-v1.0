#ifndef FRAMERT_H_
#define FRAMERT_H_

#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <vector>

#include "types.h"
#include "matrix3D.h"

class frameRT{

public:
	frameRT();
	frameRT(matrix3D,vectorOfDoubles);
	virtual ~frameRT();
			
	void setTranslation(vectorOfDoubles);
	void setRotation(matrix3D);
	
	vectorOfDoubles getTranslation();
	matrix3D getRotation();
	
	vectorOfDoubles transform(vectorOfDoubles);
	
	vectorOfDoubles sub(vectorOfDoubles,vectorOfDoubles);
	vectorOfDoubles add(vectorOfDoubles,vectorOfDoubles);
	double dot(vectorOfDoubles,vectorOfDoubles);
	vectorOfDoubles cross3(vectorOfDoubles,vectorOfDoubles);
	double norm(vectorOfDoubles);
	
private:
	vectorOfDoubles rAtom1;
	vectorOfDoubles rAtom2;
	vectorOfDoubles rAtom3;
	
	vectorOfDoubles translation;
	matrix3D rotation;
	
	
};

#endif
