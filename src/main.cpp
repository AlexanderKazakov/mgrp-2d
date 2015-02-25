#include <cstdlib>
#include <iostream>
#include <string>
#include <unistd.h>

#include <tinyxml.h>

#include "util.hpp"
#include "Stratum.hpp"
#include "Fracture.hpp"

/**
 *         Multistage Hydraulic Fracturing - 2D
 * 
 * A project for calculation of fractures' propagation 
 * with multistage hydraulic fracturing in stratum 
 * using two-dimensional approach in quasi-static linear elastic theory.
 * 
 * Used numerical method is displacement discontinuity method 
 * by Crouch (Crouch, S.L. and Starfield, A.M. (1983).
 * Boundary Element Methods in Solid Mechanics) with a few
 * improvements.
 *  
 */


/**
 * Function to load input data and parameters from xml-file.
 * @param stratum an instance of Stratum to set parameters
 * @param taskfile a file to load from
 */
void loadTask(Stratum &stratum, const char* taskfile) {
	info("Loading task from", taskfile, "...");
	TiXmlDocument *xml_file = new TiXmlDocument(taskfile);
	if (!xml_file->LoadFile()) {
		throw ("Specified taskfile is invalid");
	}
	TiXmlElement *xml_task = xml_file->FirstChildElement("task");

	TiXmlElement *xml_system = xml_task->FirstChildElement("system");
	std::string sequence = xml_system->Attribute("sequence_of_calculation");
	std::string tip = xml_system->Attribute("tip");
	std::string rotation = xml_system->Attribute("rotation");
	stratum.setSequence(sequence);
	stratum.setRotation(rotation);

	TiXmlElement *xml_stratum = xml_task->FirstChildElement("stratum");
	TiXmlElement *xml_elastic_modules = xml_stratum->FirstChildElement("elastic_modules");
	double G = atof(xml_elastic_modules->Attribute("G"));
	double nu = atof(xml_elastic_modules->Attribute("nu"));
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
		double volume = atof(xml_fracture->Attribute("volume"));
		
		TiXmlElement *initial = xml_fracture->FirstChildElement("initial");
		double x = atof(initial->Attribute("x"));
		double y = atof(initial->Attribute("y"));
		double beta = M_PI * atof(initial->Attribute("angle")) / 180;

		TiXmlElement *elements = xml_fracture->FirstChildElement("elements");
		double halfLengthOfElements = atof(elements->Attribute("half-length"));

		TiXmlElement *xml_pressure = xml_fracture->FirstChildElement("pressure");
		double a = atof(xml_pressure->Attribute("a"));
		double b = atof(xml_pressure->Attribute("b"));
		double c = atof(xml_pressure->Attribute("c"));
		std::string pressureType = (xml_pressure->Attribute("type"));

		stratum.addFracture(number, volume, x, y, beta, halfLengthOfElements,
				a, b, c, pressureType, tip);
		xml_fracture = xml_fracture->NextSiblingElement("fracture");
	}
	info("Task from", taskfile, "is loaded");
}


int main(int argc, char** argv) {
	int opt = 0;
	char *taskfile;
	bool drawDisplacement = false;
	while ((opt = getopt(argc, argv, "t:d:")) != -1) {
		switch (opt) {
			case 't': taskfile = optarg;
				break;
			case 'd': drawDisplacement = true;
				break;
			case '?': {
				info("Usage:\t./mhf-2d -t name_of_taskfile");
				exit(-1);
			}
		}
	}
	
	info("MHF-2D is running ...");
	Stratum stratum;
	try {
		loadTask(stratum, taskfile);
	} catch (const char *str) {
		print("Loading task is failed:\n", str);
		exit(-1);
	}
	stratum.calculateTask();
	stratum.visualize();
	if ( drawDisplacement )
		stratum.drawDisplacements();
	info("MHF-2D finished the work.\n");
	return 0;
}