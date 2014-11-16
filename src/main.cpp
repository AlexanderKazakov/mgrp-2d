#include <cstdlib>
#include <iostream>
#include <iostream>
#include <string>
#include <tinyxml.h>
#include <unistd.h>

#include"Stratum.hpp"
#include"Fracture.hpp"

void loadTask(Stratum &stratum, const char* taskfile);

int main(int argc, char** argv) {
	int opt = 0;
	char *taskfile;
	while ((opt = getopt(argc, argv, "h:t:")) != -1) {
		switch (opt) {
			case 'h': std::cout << "" << std::endl;
				break;
			case 't': taskfile = optarg;
				break;
			case '?': exit(-1);
		}
	}
	
	Stratum stratum;
	
	loadTask(stratum, taskfile);

	while (stratum.calculateNextFracture());

	stratum.visualize();
	
	return 0;
}

void loadTask(Stratum &stratum, const char* taskfile) {

	TiXmlDocument *xml_file = new TiXmlDocument(taskfile);
	if(!xml_file->LoadFile()) {
		std::cout << "Bad task file" << std::endl;
		return;
	}
	TiXmlElement *xml_task = xml_file->FirstChildElement("task");
	
	TiXmlElement *xml_stratum = xml_task->FirstChildElement("stratum");
	double G = atof( xml_stratum->Attribute("G") );
	double nu = atof( xml_stratum->Attribute("nu") );
	stratum.setRheology(G, nu);
	TiXmlElement *xml_ranges = xml_task->FirstChildElement("ranges");
	double Xmin = atof(xml_ranges->Attribute("Xmin"));
	double Xmax = atof(xml_ranges->Attribute("Xmax"));
	double Ymin = atof(xml_ranges->Attribute("Ymin"));
	double Ymax = atof(xml_ranges->Attribute("Ymax"));
	stratum.setRanges(Xmin, Xmax, Ymin, Ymax);
	
	TiXmlElement *xml_fractures = xml_task->FirstChildElement("fractures");
	TiXmlElement *xml_fracture = xml_fractures->FirstChildElement("fracture");
	while (xml_fracture != NULL) {
		
		int number = atoi(xml_fracture->Attribute("number"));
		TiXmlElement *initial = xml_fracture->FirstChildElement("initial");
		double x = atof(initial->Attribute("x"));
		double y = atof(initial->Attribute("y"));
		double beta = M_PI * atof(initial->Attribute("angle")) / 180;
		TiXmlElement *elements = xml_fracture->FirstChildElement("elements");
		int numOfElems = atoi(elements->Attribute("number_of_elements"));
		if (numOfElems % 2 == 0) numOfElems += 1;
		double half_length = atof(elements->Attribute("half-length"));
		TiXmlElement *xml_pressure = xml_fracture->FirstChildElement("pressure");
		double a = atof(xml_pressure->Attribute("a"));
		double b = atof(xml_pressure->Attribute("b"));
		double c = atof(xml_pressure->Attribute("c"));
		std::string pressureType = (xml_pressure->Attribute("type"));
		
		stratum.addFracture( Fracture(number, x, y, beta, half_length, numOfElems,
												a, b, c, G, nu, pressureType) );
		xml_fracture = xml_fracture->NextSiblingElement("fracture");
	}
	stratum.sortFractures();
}