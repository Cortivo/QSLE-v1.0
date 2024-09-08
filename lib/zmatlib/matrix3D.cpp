/*
 *  matrix3D.cpp
 *  gc
 *
 *  Created by Mirco Zerbetto on 7/1/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "matrix3D.h"

/***************/
/* Constructor */
/***************/

matrix3D::matrix3D()
{
	xx = 0.0;
	xy = 0.0;
	xz = 0.0;
	yx = 0.0;
	yy = 0.0;
	yz = 0.0;
	zx = 0.0;
	zy = 0.0;
	zz = 0.0;
}
matrix3D::matrix3D(double mxx, double mxy, double mxz, double myx, double myy, double myz, double mzx, double mzy, double mzz)
{
	xx = mxx;
	xy = mxy;
	xz = mxz;
	yx = myx;
	yy = myy;
	yz = myz;
	zx = mzx;
	zy = mzy;
	zz = mzz;
}


/****************/
/* Decontructor */
/****************/
matrix3D::~matrix3D()
{
}

/*********************/
/* Set entire row(s) */
/*********************/

void matrix3D::setRow1(vectorOfDoubles r)
{
	if (r.size() == 3)
	{
		xx = r[0];
		xy = r[1];
		xz = r[2];
	}
	else
	{
		std::cout << "ERROR >> matrix3D::setRow1(vectorOfDoubles) : vector size " << r.size() << " does not match 3" << std::endl;
		exit(1);
	}
	return;
}

void matrix3D::setRow2(vectorOfDoubles r)
{
	if (r.size() == 3)
	{
		yx = r[0];
		yy = r[1];
		yz = r[2];
	}
	else
	{
		std::cout << "ERROR >> matrix3D::setRow2(vectorOfDoubles) : vector size " << r.size() << " does not match 3" << std::endl;
		exit(1);
	}
	return;
}

void matrix3D::setRow3(vectorOfDoubles r)
{
	if (r.size() == 3)
	{
		zx = r[0];
		zy = r[1];
		zz = r[2];
	}
	else
	{
		std::cout << "ERROR >> matrix3D::setRow3(vectorOfDoubles) : vector size " << r.size() << " does not match 3" << std::endl;
		exit(1);
	}
	return;
}

vectorOfDoubles matrix3D::getRow1()
{
	vectorOfDoubles r;
	r.push_back(xx);
	r.push_back(xy);
	r.push_back(xz);
	return r;
}

vectorOfDoubles matrix3D::getRow2()
{
	vectorOfDoubles r;
	r.push_back(yx);
	r.push_back(yy);
	r.push_back(yz);
	return r;
}

vectorOfDoubles matrix3D::getRow3()
{
	vectorOfDoubles r;
	r.push_back(zx);
	r.push_back(zy);
	r.push_back(zz);
	return r;
}

void matrix3D::setRows(vectorOfDoubles r1, vectorOfDoubles r2, vectorOfDoubles r3)
{
	setRow1(r1);
	setRow2(r2);
	setRow3(r3);
	return;	
}

/*********************/
/* Set entire col(s) */
/*********************/

void matrix3D::setCol1(vectorOfDoubles c)
{
	if (c.size() == 3)
	{
		xx = c[0];
		yx = c[1];
		zx = c[2];
	}
	else
	{
		std::cout << "ERROR >> matrix3D::setCol1(vectorOfDoubles) : vector size " << c.size() << " does not match 3" << std::endl;
		exit(1);
	}
	return;
}

void matrix3D::setCol2(vectorOfDoubles c)
{
	if (c.size() == 3)
	{
		xy = c[0];
		yy = c[1];
		zy = c[2];
	}
	else
	{
		std::cout << "ERROR >> matrix3D::setCol2(vectorOfDoubles) : vector size " << c.size() << " does not match 3" << std::endl;
		exit(1);
	}
	return;
}

void matrix3D::setCol3(vectorOfDoubles c)
{
	if (c.size() == 3)
	{
		xz = c[0];
		yz = c[1];
		zz = c[2];
	}
	else
	{
		std::cout << "ERROR >> matrix3D::setCol3(vectorOfDoubles) : vector size " << c.size() << " does not match 3" << std::endl;
		exit(1);
	}
	return;
}

vectorOfDoubles matrix3D::getCol1()
{
	vectorOfDoubles c;
	c.push_back(xx);
	c.push_back(yx);
	c.push_back(zx);
	return c;
}

vectorOfDoubles matrix3D::getCol2()
{
	vectorOfDoubles c;
	c.push_back(xy);
	c.push_back(yy);
	c.push_back(zy);
	return c;
}

vectorOfDoubles matrix3D::getCol3()
{
	vectorOfDoubles c;
	c.push_back(xz);
	c.push_back(yz);
	c.push_back(zz);
	return c;
}

void matrix3D::setCols(vectorOfDoubles c1, vectorOfDoubles c2, vectorOfDoubles c3)
{
	setCol1(c1);
	setCol2(c2);
	setCol3(c3);
	return;	
}

/********************/
/* Transpose matrix */
/********************/

void matrix3D::transpose()
{
	vectorOfDoubles r1 = getRow1();
	vectorOfDoubles r2 = getRow2();
	vectorOfDoubles r3 = getRow3();
	setCol1(r1);
	setCol2(r2);
	setCol3(r3);
}

/********************************/
/* Matrix vector multiplication */
/********************************/

vectorOfDoubles matrix3D::multiply(vectorOfDoubles v)
{
	if (v.size() !=3)
	{
		std::cout << "ERROR >> matrix3D::multiply(vectorOfDoubles) : size of vector, " << v.size() << " does not match 3" << std::endl;
		exit(1);
	}
	vectorOfDoubles result = vectorOfDoubles(3,0.0);;
	result[0] = xx*v[0] + xy*v[1] + xz*v[2];
	result[1] = yx*v[0] + yy*v[1] + yz*v[2];
	result[2] = zx*v[0] + zy*v[1] + zz*v[2];
	return result;
}

/********************************/
/* matrix matrix multiplication */
/********************************/

void matrix3D::multiply(matrix3D a, matrix3D b)
{
	xx = a.xx*b.xx + a.xy*b.yx + a.xz*b.zx;
	xy = a.xx*b.xy + a.xy*b.yy + a.xz*b.zy;
	xz = a.xx*b.xz + a.xy*b.yz + a.xz*b.zz;

	yx = a.yx*b.xx + a.yy*b.yx + a.yz*b.zx;
	yy = a.yx*b.xy + a.yy*b.yy + a.yz*b.zy;
	yz = a.yx*b.xz + a.yy*b.yz + a.yz*b.zz;

	zx = a.zx*b.xx + a.zy*b.yx + a.zz*b.zx;
	zy = a.zx*b.xy + a.zy*b.yy + a.zz*b.zy;
	zz = a.zx*b.xz + a.zy*b.yz + a.zz*b.zz;

	return;
}
