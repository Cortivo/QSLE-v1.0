/*
 *  config.h
 *  gc
 *
 *  Created by Mirco Zerbetto on 6/30/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#include <cstdio>
#include <cstdlib>
#include <string>

class config{
	
public:
	
	config();
	virtual ~config();
	
	void setSesto(std::string);
	std::string getSesto();
	
	void setVdwFile(std::string);
	std::string getVdwFile();
	
	void setAminoAcidsFile(std::string);
	std::string getAminoAcidsFile();
		
	void setMassFile(std::string);
	std::string getMassFile();

private:
	
	std::string sesto;
	std::string vdwFile;
	std::string aminoAcidsFile;
	std::string massFile;

};

#endif
