/*
 *  config.cpp
 *  gc
 *
 *  Created by Mirco Zerbetto on 6/30/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "config.h"

/***************/
/* Constructor */
/***************/

config::config()
{	
}

/*****************/
/* Deconstructor */
/*****************/

config::~config()
{
}

/*****************/
/* SeStO program */
/*****************/

void config::setSesto(std::string p)
{
	sesto = p;
	return;
}
std::string config::getSesto()
{
	return sesto;
}

/****************/
/* vdw.dat file */
/****************/

void config::setVdwFile(std::string f)
{
	vdwFile = f;
	return;
}

std::string config::getVdwFile()
{
	return vdwFile;
}

/***********************/
/* Aminoacids database */
/***********************/

void config::setAminoAcidsFile(std::string f)
{
	aminoAcidsFile = f;
	return;
}

std::string config::getAminoAcidsFile()
{
	return aminoAcidsFile;
}

/*******************/
/* Masses database */
/*******************/

void config::setMassFile(std::string f)
{
	massFile = f;
	return;
}

std::string config::getMassFile()
{
	return massFile;
}
