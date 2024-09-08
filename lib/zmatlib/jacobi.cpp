/*
 *  jacobi.cpp
 *  gc
 *
 *  Created by Mirco Zerbetto on 7/24/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "jacobi.h"

/****************/
/* Constructors */
/****************/

jacobi::jacobi()
{
	haveMatrix = false;
	a.clear();
	v.clear();
	d.clear();
	nrot = 0;
	maxStep = 50;
}

jacobi::jacobi(vectorOfDoubles m)
{
	setMatrix(m);
	maxStep = 50;
}

/**************/
/* Destructor */
/**************/

jacobi::~jacobi()
{
}

/*********************************/
/* Pass the matrix to the object */
/*********************************/

void jacobi::setMatrix(vectorOfDoubles m)
{
	if (m.size() > 0)
	{
		haveMatrix = true;
		a = m;
		v = vectorOfDoubles(a.size(),0.0);
		d = vectorOfDoubles((int)sqrt(a.size()),0.0);
		nrot = 0;
	}
	else
	{
		haveMatrix = false;
		a.clear();
		v.clear();
		d.clear();
		nrot = 0;		
		std::cout << "WARNING >> jacobi::setMatrix(vectorOfDoubles) : passed matrix is not allocated. Matrix not present in jacobi() object" << std::endl;
	}	
}

void jacobi::setMatrix(matrix3D m)
{
	vectorOfDoubles elements;
	elements.push_back(m.xx);
	elements.push_back(m.xy);
	elements.push_back(m.xz);
	elements.push_back(m.yx);
	elements.push_back(m.yy);
	elements.push_back(m.yz);
	elements.push_back(m.zx);
	elements.push_back(m.zy);
	elements.push_back(m.zz);
	setMatrix(elements);
	return;	
}

/*******************************/
/* Overwrite default max steps */
/*******************************/

void jacobi::setMaxStep(int n)
{
	maxStep = n;
	return;
}

/******************/
/* Return results */
/******************/

vectorOfDoubles jacobi::getMatrix()
{
	if (!haveMatrix)
	{
		std::cout << "ERROR >> jacobi::gatMatrix() : matrix was not allocated"<<std::endl;
		exit(1);
	}

	return a;
}

vectorOfDoubles jacobi::getEigenValues()
{
	if (!haveMatrix)
	{
		std::cout << "ERROR >> jacobi::gatEigenValues() : matrix was not allocated"<<std::endl;
		exit(1);
	}

	return d;
}

vectorOfDoubles jacobi::getEigenVectors()
{
	if (!haveMatrix)
	{
		std::cout << "ERROR >> jacobi::gatEigenVectors() : matrix was not allocated"<<std::endl;
		exit(1);
	}
	
	return v;
}

matrix3D jacobi::getMatrix3D()
{
	if (!haveMatrix)
	{
		std::cout << "ERROR >> jacobi::getMatrix3D() : matrix was not allocated"<<std::endl;
		exit(1);
	}
	
	if (a.size()==9)
	{
		matrix3D m;
		m.xx = a[0]; m.xy = a[1]; m.xz = a[2];
		m.yx = a[3]; m.yy = a[4]; m.yz = a[5];
		m.zx = a[6]; m.zy = a[7]; m.zz = a[8];
		return m;
	}
	else
	{
		std::cout << "ERROR: jacobi::getMatrix3D() : cannot return a matrix3D object for n != 3 (present value is " << d.size() << ")" << std::endl;
		exit(1);
	}
}

matrix3D jacobi::getEigenVectors3D()
{
	if (!haveMatrix)
	{
		std::cout << "ERROR >> jacobi::getEigenVectors3D() : matrix was not allocated"<<std::endl;
		exit(1);
		}
		
		if (a.size()==9)
		{
			matrix3D m;
			m.xx = v[0]; m.xy = v[1]; m.xz = v[2];
			m.yx = v[3]; m.yy = v[4]; m.yz = v[5];
			m.zx = v[6]; m.zy = v[7]; m.zz = v[8];
			return m;
		}
		else
		{
			std::cout << "ERROR: jacobi::getEigenVectors3D() : cannot return a matrix3D object for n != 3 (present value is " << d.size() << ")" << std::endl;
			exit(1);
		}
}
		
matrix3D jacobi::getEigenValues3D()
{
	if (!haveMatrix)
	{
		std::cout << "ERROR >> jacobi::getEigenValues3D() : matrix was not allocated"<<std::endl;
		exit(1);
	}
	
	if (a.size()==9)
	{
		matrix3D m;
		m.xx = d[0]; m.xy = 0.0 ; m.xz = 0.0 ;
		m.yx = 0.0 ; m.yy = d[1]; m.yz = 0.0 ;
		m.zx = 0.0 ; m.zy = 0.0 ; m.zz = d[2];
		return m;
	}
	else
	{
		std::cout << "ERROR: jacobi::getEigenValues3D() : cannot return a matrix3D object for n != 3 (present value is " << d.size() << ")" << std::endl;
		exit(1);
	}
}

int jacobi::getNrot()
{
	return nrot;
}

/*******************/
/* Diagonalization */
/*******************/

/////////////////////////////////////////////////////////
// From: NUMERICAL RECIPES IN C, chap. 11, pp. 463-469 //
/////////////////////////////////////////////////////////

#define TOL 1.0e-15

void jacobi::diagonalize()
{
	if (!haveMatrix)
	{
		std::cout << "ERROR >> jacobi::diagonalize() : matrix was not allocated"<<std::endl;
		exit(1);
	}
	
	int n = d.size();
	int j, iq, ip, i;
	double tresh, theta, tau, t, sm, s, h, g, c;

	vectorOfDoubles b(n,0.0);
	vectorOfDoubles z(n,0.0);	
	for (ip = 0; ip < n; ip++)
	{
		// set v to the identity matrix
		for (iq = 0; iq < n; iq++)
			v[ip*n+iq] = 0.0;
		v[ip*n+ip] = 1.0;
		// copy the diagonal of a into b and d
		b[ip] = d[ip] = a[ip*n+ip];
	}

	nrot = 0;
	for (i = 0; i < maxStep; i++)
	{
		// sum off-diagonal elements
		sm = 0.0;
		for (ip = 0; ip < n-1; ip++)
		{
			for (iq = ip+1; iq < n; iq++)
				sm += fabs(a[ip*n+iq]);
		}
		
		// is the matrix diagonal?
		if (fabs(sm) < TOL)
		{
			// Check if the rotation is proper and
			// correct the matrix if not

  			double detV = v[0] * (v[4] * v[8] - v[5] * v[7]) - v[1] * (v[3] * v[8] - v[5] * v[6]) + v[2] * (v[3] * v[7] - v[4] * v[6]);
  
    			//std::cout << "++++   determinant is now: " << detV << std::endl;
			if (detV < 0.0)
			{
				std::cout << std::endl << "Jacobi diagonalization: eigenvectors matrix has determinant -1. Restoring the right-handed rotation changing sign to Z" << std::endl;
				v[0] = -v[0];
				v[1] = -v[1];
				v[2] = -v[2];
				v[3] = -v[3];
				v[4] = -v[4];
				v[5] = -v[5];
				v[6] = -v[6];
				v[7] = -v[7];
				v[8] = -v[8];
  				detV = v[0] * (v[4] * v[8] - v[5] * v[7]) - v[1] * (v[3] * v[8] - v[5] * v[6]) + v[2] * (v[3] * v[7] - v[4] * v[6]);
    				std::cout << "   determinant is now: " << detV << std::endl;
			}
			return;
		}
		
		// calculate tresh
		
		if (i<4)
			tresh = 0.2 * sm / (double)(n*n);
		else
			tresh = 0.0;
		
		// perform rotation
		
		for (ip = 0; ip < n-1; ip++)
		{
			for (iq = ip+1; iq < n; iq++)
			{
				g = 100.0*fabs(a[ip*n+iq]);
				
				if (i>4 && (fabs(d[ip])+g)==fabs(d[ip]) && (fabs(d[iq])+g)==fabs(d[iq]))
					a[ip*n+iq] = 0.0;
				else if (fabs(a[ip*n+iq]) > tresh)
				{
					h = d[iq] - d[ip];
					if( (fabs(h)+g) == fabs(h) )
						t = a[ip*n+iq] / h;
					else
					{
						theta = 0.5 *h / a[ip*n+iq];
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta));
						if (theta < 0.0) t = -t;
					}
					
					c = 1.0 / sqrt(1.0 + t*t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[ip*n+iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip*n+iq] = 0.0;
					
					for (j = 0; j < ip; j++)
						ROTATE(0,j,ip,j,iq,n,s,tau);
					for (j = ip+1; j < iq; j++)
						ROTATE(0,ip,j,j,iq,n,s,tau);
					for (j = iq+1; j < n; j++)
						ROTATE(0,ip,j,iq,j,n,s,tau);
					for (j = 0; j < n; j++)
						ROTATE(1,j,ip,j,iq,n,s,tau);
					
					nrot++;
						
				} // else if
				
			}// iq
		}// ip
		
		for (ip = 0; ip < n; ip++)
		{
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
		
	}// i
	
	std::cout  << "WARNING >> jacobi::diagonalize() : Max n iterations (nrot = " << maxStep << ") reached" << std::endl;
	return;
}

void jacobi::ROTATE(int mat, int i, int j, int k, int l, int n, double s, double tau)
{
	switch(mat)
	{
		case 0:
		{
			double g = a[i*n+j];
			double h = a[k*n+l];
			a[i*n+j] = g - s*(h+g*tau);
			a[k*n+l] = h + s*(g-h*tau);
			break;
		}
		default:
		{
			double g = v[i*n+j];
			double h = v[k*n+l];
			v[i*n+j] = g - s*(h+g*tau);
			v[k*n+l] = h + s*(g-h*tau);
			break;
		}
	}
	return;	
}

/************************/
/* Print debugging info */
/************************/

void jacobi::debug()
{
	std::cout << std::endl << "**** MATRIX ****" << std::endl;
	for (unsigned int i = 0; i < d.size(); i++)
	{
		for (unsigned int j = 0; j < d.size(); j++)
			std::cout << a[i*d.size()+j] << "\t";
		std::cout << std::endl;
	}

	std::cout << std::endl << "**** EIGENVECTORS ****" << std::endl;
	for (unsigned int i = 0; i < d.size(); i++)
	{
		for (unsigned int j = 0; j < d.size(); j++)
			std::cout << v[i*d.size()+j] << "\t";
		std::cout << std::endl;
	}
	
	std::cout << std::endl << "**** EIGENVALUES ****" << std::endl;

	for (unsigned int i = 0; i < d.size(); i++)
		std::cout << d[i] << std::endl;

	std::cout << std::endl << "**** NROT ****" << std::endl;
	std::cout << nrot << " / " << maxStep << std::endl << std::endl;
	
	return;
}

void jacobi::reorder()
{
	int i,j,k,n;
	double p;
	n = d.size();
	for (i = 0; i < n-1; i++)
	{
		p = d[k=i];
		for (j = i+1; j < n; j++)
			if (d[j] <= p) p = d[k=j];
		if (k != i)
		{
			d[k] = d[i];
			d[i] = p;
			for (j = 0; j < n; j++)
			{
				p = v[j*n+i];
				v[j*n+i] = v[j*n+k];
				v[j*n+k] = p;
			}
		}
	}
}

