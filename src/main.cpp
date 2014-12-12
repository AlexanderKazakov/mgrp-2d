#include <cstdlib>
#include <iostream>
#include <iostream>
#include <string>
#include <tinyxml.h>
#include <unistd.h>

#include"Stratum.hpp"
#include"Fracture.hpp"

/**
 * A project for calculation of fractures' propagation 
 * with multistage hydraulic fracturing in stratum 
 * using 2D-approach in quasi-static linear elastic theory.
 * 
 * Used numerical method is displacement discontinuity method 
 * by Crouch (Crouch, S.L. and Starfield, A.M. (1983).
 * Boundary Element Methods in Solid Mechanics) with a few
 * improvements.
 *  
 */


void loadTask(Stratum &stratum, const char* taskfile);

int main(int argc, char** argv) {
	int opt = 0;
	char *taskfile;
	while ((opt = getopt(argc, argv, "t:")) != -1) {
		switch (opt) {
			case 't': taskfile = optarg;
				break;
			case '?': {
				std::cout << "Usage:	./mgrp-2d -t name_of_taskfile\n";
				exit(-1);
			}
		}
	}
	
	Stratum stratum;
	loadTask(stratum, taskfile);
	stratum.calculate();
	stratum.visualize();
	
	return 0;
}

void loadTask(Stratum &stratum, const char* taskfile) {

	TiXmlDocument *xml_file = new TiXmlDocument(taskfile);
	if(!xml_file->LoadFile()) {
		std::cout << "Usage:	./mgrp-2d -t name_of_taskfile\n" << 
									"Specified taskfile is invalid\n";
		exit(-1);
	}
	TiXmlElement *xml_task = xml_file->FirstChildElement("task");
	
	TiXmlElement *xml_system = xml_task->FirstChildElement("system");
	std::string tip = xml_system->Attribute("tip");
	std::string rotation = xml_system->Attribute("rotation");
	
	TiXmlElement *xml_stratum = xml_task->FirstChildElement("stratum");
	TiXmlElement *xml_elastic_modules = xml_stratum->FirstChildElement("elastic_modules");
	double G = atof( xml_elastic_modules->Attribute("G") );
	double nu = atof( xml_elastic_modules->Attribute("nu") );
	stratum.setRheology(G, nu);

	TiXmlElement *xml_external_stresses = xml_stratum->FirstChildElement("external_stresses");
	double Sxx = atof(xml_external_stresses->Attribute("Sxx"));
	double Sxy = atof(xml_external_stresses->Attribute("Sxy"));
	double Syy = atof(xml_external_stresses->Attribute("Syy"));
	stratum.setStresses(Sxx, Sxy, Syy);

	TiXmlElement *xml_ranges = xml_task->FirstChildElement("ranges");
	double Xmin = atof(xml_ranges->Attribute("Xmin"));
	double Xmax = atof(xml_ranges->Attribute("Xmax"));
	double Ymin = atof(xml_ranges->Attribute("Ymin"));
	double Ymax = atof(xml_ranges->Attribute("Ymax"));
	stratum.setRanges(Xmin, Xmax, Ymin, Ymax);
	
	TiXmlElement *xml_fractures = xml_task->FirstChildElement("fractures");
	TiXmlElement *xml_fracture = xml_fractures->FirstChildElement("fracture");
	int numberOfFractures = 0;
	while (xml_fracture != NULL) {
		numberOfFractures++;
		xml_fracture = xml_fracture->NextSiblingElement("fracture");
	}
	stratum.reserve(numberOfFractures);
	xml_fracture = xml_fractures->FirstChildElement("fracture");
	while (xml_fracture != NULL) {
		int number = atoi(xml_fracture->Attribute("number"));
		
		TiXmlElement *initial = xml_fracture->FirstChildElement("initial");
		double x = atof(initial->Attribute("x"));
		double y = atof(initial->Attribute("y"));
		double beta = M_PI * atof(initial->Attribute("angle")) / 180;
		
		TiXmlElement *elements = xml_fracture->FirstChildElement("elements");
		int numOfElems = atoi(elements->Attribute("number_of_elements"));
		// Number of elements in fracture is always odd because 
		// it doesn't matter but is very helpful
		if (numOfElems % 2 == 0) numOfElems += 1;
		double half_length = atof(elements->Attribute("half-length"));
		
		TiXmlElement *xml_pressure = xml_fracture->FirstChildElement("pressure");
		double a = atof(xml_pressure->Attribute("a"));
		double b = atof(xml_pressure->Attribute("b"));
		double c = atof(xml_pressure->Attribute("c"));
		std::string pressureType = (xml_pressure->Attribute("type"));

		stratum.addFracture(number, x, y, beta, half_length,
						numOfElems, a, b, c, pressureType, tip, rotation);
		xml_fracture = xml_fracture->NextSiblingElement("fracture");

	}
}
