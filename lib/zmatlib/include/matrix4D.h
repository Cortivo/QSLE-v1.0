#ifndef MATRIX4D
#define MATRIX4D

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <iomanip>

#include "types.h"

class matrix4D{

public:
	matrix4D();
	matrix4D(double,double,double,double,
		 double,double,double,double,
		 double,double,double,double,
		 double,double,double,double);
	virtual ~matrix4D();

	void print(void);
	
	void setRow(vectorOfDoubles, int);
	
	vectorOfDoubles getRow(int);

	void setCol(vectorOfDoubles, int);
	
	vectorOfDoubles getCol(int);
	
	void setRows(vectorOfDoubles,vectorOfDoubles,vectorOfDoubles,vectorOfDoubles);
	void setCols(vectorOfDoubles,vectorOfDoubles,vectorOfDoubles,vectorOfDoubles);

	void transpose();
	vectorOfDoubles multiply(vectorOfDoubles);
	void multiply(matrix4D, matrix4D);
	
	double M11, M12, M13, M14;
	double M21, M22, M23, M24;
	double M31, M32, M33, M34;
	double M41, M42, M43, M44;
	
private:
	
};

#endif
