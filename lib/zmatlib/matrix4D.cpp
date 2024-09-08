/*
 *  matrix4D.cpp
 *  gc
 *
 *  Created by Mirco Zerbetto on 7/1/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "matrix4D.h"

/***************/
/* Constructor */
/***************/

matrix4D::matrix4D()
{
	M11 = 0.0;
	M12 = 0.0;
	M13 = 0.0;
	M14 = 0.0;
	M21 = 0.0;
	M22 = 0.0;
	M23 = 0.0;
	M24 = 0.0;
	M31 = 0.0;
	M32 = 0.0;
	M33 = 0.0;
	M34 = 0.0;
	M41 = 0.0;
	M42 = 0.0;
	M43 = 0.0;
	M44 = 0.0;
}
matrix4D::matrix4D(double a11, double a12, double a13, double a14,
		   double a21, double a22, double a23, double a24,
		   double a31, double a32, double a33, double a34,
		   double a41, double a42, double a43, double a44)
{
	M11 = a11;
	M12 = a12;
	M13 = a13;
	M14 = a14;
	M21 = a21;
	M22 = a22;
	M23 = a23;
	M24 = a24;
	M31 = a31;
	M32 = a32;
	M33 = a33;
	M34 = a34;
	M41 = a41;
	M42 = a42;
	M43 = a43;
	M44 = a44;
}


/****************/
/* Decontructor */
/****************/
matrix4D::~matrix4D()
{
}

/****************/
/* Print matrix */
/****************/

void matrix4D::print(void)
{
	std::cout << std::endl;
	std::cout << M11 << "\t" << M12 << "\t" << M13 << "\t" << M14 << std::endl;
	std::cout << M21 << "\t" << M22 << "\t" << M23 << "\t" << M24 << std::endl;
	std::cout << M31 << "\t" << M32 << "\t" << M33 << "\t" << M34 << std::endl;
	std::cout << M41 << "\t" << M42 << "\t" << M43 << "\t" << M44 << std::endl;
	std::cout << std::endl;
	return;
}

/*********************/
/* Set entire row(s) */
/*********************/

void matrix4D::setRow(vectorOfDoubles r, int idx)
{
	if (r.size() != 4)
	{
		std::cout << "ERROR >> matrix4D::setRow(vectorOfDoubles, int) : vector size " << r.size() << " does not match 4" << std::endl;
		exit(1);
	}
	switch (idx)
	{
		case 1:
		{
			M11 = r[0];
			M12 = r[1];
			M13 = r[2];
			M14 = r[3];
		}
		case 2:
		{
			M21 = r[0];
			M22 = r[1];
			M23 = r[2];
			M24 = r[3];
		}
		case 3:
		{
			M31 = r[0];
			M32 = r[1];
			M33 = r[2];
			M34 = r[3];
		}
		case 4:
		{
			M41 = r[0];
			M42 = r[1];
			M43 = r[2];
			M44 = r[3];
		}
		default:
		{
			std::cout << "ERROR >> matrix4D::setRow(vectorOfDoubles, int) : row number " << idx << " must be in [1,4]" << std::endl;
			exit(1);
		}
	}
	return;
}

vectorOfDoubles matrix4D::getRow(int idx)
{
	vectorOfDoubles r;
	switch (idx)
	{
		case 1:
		{
			r.push_back(M11);
			r.push_back(M12);
			r.push_back(M13);
			r.push_back(M14);
		}
		case 2:
		{
			r.push_back(M21);
			r.push_back(M22);
			r.push_back(M23);
			r.push_back(M24);
		}
		case 3:
		{
			r.push_back(M31);
			r.push_back(M32);
			r.push_back(M33);
			r.push_back(M34);
		}
		case 4:
		{
			r.push_back(M41);
			r.push_back(M42);
			r.push_back(M43);
			r.push_back(M44);
		}
		default:
		{
			std::cout << "ERROR >> matrix4D::setRow(vectorOfDoubles, int) : row number " << idx << " must be in [1,4]" << std::endl;
			exit(1);
		}
	}
	return r;
}

void matrix4D::setRows(vectorOfDoubles r1, vectorOfDoubles r2, vectorOfDoubles r3, vectorOfDoubles r4)
{
	setRow(r1, 1);
	setRow(r2, 2);
	setRow(r3, 3);
	setRow(r4, 4);
	return;	
}

/*********************/
/* Set entire col(s) */
/*********************/

void matrix4D::setCol(vectorOfDoubles c, int idx)
{
	if (c.size() != 4)
	{
		std::cout << "ERROR >> matrix4D::setCol(vectorOfDoubles, int) : vector size " << c.size() << " does not match 4" << std::endl;
		exit(1);
	}
	switch (idx)
	{
		case 1:
		{
			M11 = c[0];
			M21 = c[1];
			M31 = c[2];
			M41 = c[3];
		}
		case 2:
		{
			M12 = c[0];
			M22 = c[1];
			M32 = c[2];
			M42 = c[3];
		}
		case 3:
		{
			M13 = c[0];
			M23 = c[1];
			M33 = c[2];
			M43 = c[3];
		}
		case 4:
		{
			M14 = c[0];
			M24 = c[1];
			M34 = c[2];
			M44 = c[3];
		}
		default:
		{
			std::cout << "ERROR >> matrix4D::setCol(vectorOfDoubles, int) : column number " << idx << " must be in [1,4]" << std::endl;
			exit(1);
		}
	}
	return;
}

vectorOfDoubles matrix4D::getCol(int idx)
{
	vectorOfDoubles c;
	switch (idx)
	{
		case 1:
		{
			c.push_back(M11);
			c.push_back(M21);
			c.push_back(M31);
			c.push_back(M41);
		}
		case 2:
		{
			c.push_back(M12);
			c.push_back(M22);
			c.push_back(M32);
			c.push_back(M42);
		}
		case 3:
		{
			c.push_back(M13);
			c.push_back(M23);
			c.push_back(M33);
			c.push_back(M43);
		}
		case 4:
		{
			c.push_back(M14);
			c.push_back(M24);
			c.push_back(M34);
			c.push_back(M44);
		}
		default:
		{
			std::cout << "ERROR >> matrix4D::setCol(vectorOfDoubles, int) : column number " << idx << " must be in [1,4]" << std::endl;
			exit(1);
		}
	}
	return c;
}

void matrix4D::setCols(vectorOfDoubles c1, vectorOfDoubles c2, vectorOfDoubles c3, vectorOfDoubles c4)
{
	setCol(c1, 1);
	setCol(c2, 2);
	setCol(c3, 3);
	setCol(c4, 4);
	return;	
}

/********************/
/* Transpose matrix */
/********************/

void matrix4D::transpose()
{
	vectorOfDoubles r1 = getRow(1);
	vectorOfDoubles r2 = getRow(2);
	vectorOfDoubles r3 = getRow(3);
	vectorOfDoubles r4 = getRow(4);
	setCol(r1, 1);
	setCol(r2, 2);
	setCol(r3, 3);
	setCol(r4, 4);
}

/********************************/
/* Matrix vector multiplication */
/********************************/

vectorOfDoubles matrix4D::multiply(vectorOfDoubles v)
{
	if (v.size() != 4)
	{
		std::cout << "ERROR >> matrix4D::multiply(vectorOfDoubles) : size of vector, " << v.size() << " does not match 4" << std::endl;
		exit(1);
	}
	vectorOfDoubles result = vectorOfDoubles(4, 0.0);
	result[0] = M11 * v[0] + M12 * v[1] + M13 * v[2] + M14 * v[3];
	result[1] = M21 * v[0] + M22 * v[1] + M23 * v[2] + M24 * v[3];
	result[2] = M31 * v[0] + M32 * v[1] + M33 * v[2] + M34 * v[3];
	result[3] = M41 * v[0] + M42 * v[1] + M43 * v[2] + M44 * v[3];
	return result;
}

/********************************/
/* matrix matrix multiplication */
/********************************/

void matrix4D::multiply(matrix4D a, matrix4D b)
{
	double R11 = a.M11 * b.M11 + a.M12 * b.M21 + a.M13 * b.M31 + a.M14 * b.M41;
	double R12 = a.M11 * b.M12 + a.M12 * b.M22 + a.M13 * b.M32 + a.M14 * b.M42;
	double R13 = a.M11 * b.M13 + a.M12 * b.M23 + a.M13 * b.M33 + a.M14 * b.M43;
	double R14 = a.M11 * b.M14 + a.M12 * b.M24 + a.M13 * b.M34 + a.M14 * b.M44;

	double R21 = a.M21 * b.M11 + a.M22 * b.M21 + a.M23 * b.M31 + a.M24 * b.M41;
	double R22 = a.M21 * b.M12 + a.M22 * b.M22 + a.M23 * b.M32 + a.M24 * b.M42;
	double R23 = a.M21 * b.M13 + a.M22 * b.M23 + a.M23 * b.M33 + a.M24 * b.M43;
	double R24 = a.M21 * b.M14 + a.M22 * b.M24 + a.M23 * b.M34 + a.M24 * b.M44;

	double R31 = a.M31 * b.M11 + a.M32 * b.M21 + a.M33 * b.M31 + a.M34 * b.M41;
	double R32 = a.M31 * b.M12 + a.M32 * b.M22 + a.M33 * b.M32 + a.M34 * b.M42;
	double R33 = a.M31 * b.M13 + a.M32 * b.M23 + a.M33 * b.M33 + a.M34 * b.M43;
	double R34 = a.M31 * b.M14 + a.M32 * b.M24 + a.M33 * b.M34 + a.M34 * b.M44;

	double R41 = a.M41 * b.M11 + a.M42 * b.M21 + a.M43 * b.M31 + a.M44 * b.M41;
	double R42 = a.M41 * b.M12 + a.M42 * b.M22 + a.M43 * b.M32 + a.M44 * b.M42;
	double R43 = a.M41 * b.M13 + a.M42 * b.M23 + a.M43 * b.M33 + a.M44 * b.M43;
	double R44 = a.M41 * b.M14 + a.M42 * b.M24 + a.M43 * b.M34 + a.M44 * b.M44;

	M11 = R11;
	M12 = R12;
	M13 = R13;
	M14 = R14;
                
	M21 = R21;
	M22 = R22;
	M23 = R23;
	M24 = R24;
                
	M31 = R31;
	M32 = R32;
	M33 = R33;
	M34 = R34;
                
	M41 = R41;
	M42 = R42;
	M43 = R43;
	M44 = R44;

	return;
}
