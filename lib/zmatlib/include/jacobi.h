/*
 *  jacobi.h
 *  gc
 *
 *  Created by Mirco Zerbetto on 7/24/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef JACOBI_H_
#define JACOBI_H_

#include <cstdlib>
#include <math.h>
#include <vector>

#include "matrix3D.h"

class jacobi
	{
	public:
		jacobi();
		jacobi(vectorOfDoubles);
		virtual ~jacobi();
		
		void setMatrix(vectorOfDoubles);
		void setMatrix(matrix3D);
		void setMaxStep(int);
		
		void diagonalize();
		void reorder();
		
		vectorOfDoubles getMatrix();
		vectorOfDoubles getEigenValues();
		vectorOfDoubles getEigenVectors();

		matrix3D getMatrix3D();
		matrix3D getEigenValues3D();
		matrix3D getEigenVectors3D();

		int getNrot();
		
		void debug();
		
	private:
		bool haveMatrix;
		int nrot, maxStep;
		vectorOfDoubles a;
		vectorOfDoubles d;
		vectorOfDoubles v;
		void ROTATE(int, int, int, int, int, int, double, double);
		
	};

#endif
