#ifndef ZMAT_TYPES_H
#define ZMAT_TYPES_H

#include <vector>
#include <string>

typedef std::vector<int> vectorOfIntegers;
typedef vectorOfIntegers::iterator intit;
typedef std::vector<std::string> vectorOfStrings;
typedef std::vector<double> vectorOfDoubles;

struct zmatline {
	int A, B, C;
	double d, theta, phi;
};

struct frag {
	int isref, next, nt, isTerminal;
	int natoms, *atoms;
};

#endif // ZMAT_TYPES
