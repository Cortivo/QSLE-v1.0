/*
 *  matrix3D.h
 *  gc
 *
 *  Created by Mirco Zerbetto on 7/1/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MATRIX3D_H_
#define MATRIX3D_H_

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <iomanip>

#include "types.h"

class matrix3D{

public:
	matrix3D();
	matrix3D(double,double,double,double,double,double,double,double,double);
	virtual ~matrix3D();
	
	void setRow1(vectorOfDoubles);
	void setRow2(vectorOfDoubles);
	void setRow3(vectorOfDoubles);
	
	vectorOfDoubles getRow1();
	vectorOfDoubles getRow2();
	vectorOfDoubles getRow3();

	void setCol1(vectorOfDoubles);
	void setCol2(vectorOfDoubles);
	void setCol3(vectorOfDoubles);
	
	vectorOfDoubles getCol1();
	vectorOfDoubles getCol2();
	vectorOfDoubles getCol3();
	
	void setRows(vectorOfDoubles,vectorOfDoubles,vectorOfDoubles);
	void setCols(vectorOfDoubles,vectorOfDoubles,vectorOfDoubles);

	void transpose();
	vectorOfDoubles multiply(vectorOfDoubles);
	void multiply(matrix3D, matrix3D);
	
	double xx, xy, xz;
	double yx, yy, yz;
	double zx, zy, zz;
	
private:
	
};

#endif
