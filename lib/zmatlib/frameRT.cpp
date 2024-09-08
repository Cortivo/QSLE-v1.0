#include "frameRT.h"

/***************/
/* Constructor */
/***************/

frameRT::frameRT()
{
	translation = vectorOfDoubles(3,0.0);
	rotation.xx = 1.0; rotation.xy = 0.0; rotation.xz = 0.0;
	rotation.yx = 0.0; rotation.yy = 1.0; rotation.yz = 0.0;
	rotation.zx = 0.0; rotation.zy = 0.0; rotation.zz = 1.0;
}

frameRT::frameRT(matrix3D R, vectorOfDoubles T)
{
	if (T.size() == 3)
		translation = T;
	else
	{
		std::cout << "ERROR >> frameRT(matrix3D, vectorOfDoubles) : size of vector must be 3, not " << T.size() << std::endl;
		exit(1);
	}
	
	rotation = R;
}

/*****************/
/* Deconstructor */
/*****************/

frameRT::~frameRT()
{
}

/************************/
/* Set roto-translation */
/************************/

void frameRT::setTranslation(vectorOfDoubles t)
{
	if (t.size() != 3)
	{
		std::cout << "ERROR >> frameRT::setTranslation(vectorOfDoubles) : translation vector must be of size 3, not " << t.size() << std::endl;
		exit(1);
	}
	translation = t;
	return;
}

void frameRT::setRotation(matrix3D r)
{
	rotation = r;
	return;
}

/************************/
/* Get roto-translation */
/************************/

vectorOfDoubles frameRT::getTranslation()
{
	return translation;
}

matrix3D frameRT::getRotation()
{
	return rotation;
}

/********************************************************/
/* Transform a vector using the stored roto-translation */
/********************************************************/

vectorOfDoubles frameRT::transform(vectorOfDoubles v)
{
	vectorOfDoubles vNew;
    vNew = sub(v,translation);
	vNew = rotation.multiply(vNew);
	

	return vNew;
}

/******************************/
/* Vector - Vector operations */
/******************************/

vectorOfDoubles frameRT::sub(vectorOfDoubles a, vectorOfDoubles b)
{
	if (a.size() != b.size())
	{
		std::cout << "ERROR >> frameRT::sub(vectorOfDoubles,vectorOfDoubles) : size of first vector, " << a.size() << ", does not match size of second vector, " << b.size() << std::endl;
		exit(1);
	}
	for (unsigned int i = 0; i < a.size(); i++)
		a.at(i) -= b.at(i);
	return a;
}

vectorOfDoubles frameRT::add(vectorOfDoubles a, vectorOfDoubles b)
{
	if (a.size() != b.size())
	{
		std::cout << "ERROR >> frameRT::add(vectorOfDoubles,vectorOfDoubles) : size of first vector, " << a.size() << ", does not match size of second vector, " << b.size() << std::endl;
		exit(1);
	}
	for (unsigned int i = 0; i < a.size(); i++)
		a.at(i) += b.at(i);
	return a;	
}

double frameRT::dot(vectorOfDoubles a,vectorOfDoubles b)
{
	if (a.size() != b.size())
	{
		std::cout << "ERROR >> frameRT::dot(vectorOfDoubles,vectorOfDoubles) : size of first vector, " << a.size() << ", does not match size of second vector, " << b.size() << std::endl;
		exit(1);
	}
	double d = 0.0;
	for (unsigned int i = 0; i < a.size(); i++)
		d += a.at(i)*b.at(i);
	return d;
}

vectorOfDoubles frameRT::cross3(vectorOfDoubles a, vectorOfDoubles b)
{
	if (a.size() != b.size())
	{
		std::cout << "ERROR >> frameRT::corss3(vectorOfDoubles,vectorOfDoubles) : size of first vector, " << a.size() << ", does not match size of second vector, " << b.size() << std::endl;
		exit(1);
	}
	if (a.size() != 3)
	{
		std::cout << "ERROR >> frameRT::corss3(vectorOfDoubles,vectorOfDoubles) : the funtcion calculates cross product of only 3D arrays" << std::endl;
		exit(1);		
	}
	vectorOfDoubles c(3,0.0);
	c.at(0) = a.at(1)*b.at(2) - a.at(2)*b.at(1);
	c.at(1) = a.at(2)*b.at(0) - a.at(0)*b.at(2);
	c.at(2) = a.at(0)*b.at(1) - a.at(1)*b.at(0);
	return c;	
}

double frameRT::norm(vectorOfDoubles v)
{
	double n = sqrt(dot(v,v));
	return n;
}
