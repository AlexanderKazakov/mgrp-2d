#include "Engine.hpp"

Engine::Engine() {
}

Engine::~Engine() {
}

void Engine::calculate() {
	while (stratum.calculateNextFracture());
}

void Engine::visualize() {
	stratum.visualize();
}

void Engine::loadTask() {

	TiXmlDocument *xml_file = new TiXmlDocument("task.xml");
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
		double pressure = atof(xml_pressure->Attribute("value"));
		std::string pressureType = (xml_pressure->Attribute("type"));
		
		stratum.addFracture( Fracture(number, x, y, beta, half_length, numOfElems,
											pressure, G, nu, pressureType) );
		xml_fracture = xml_fracture->NextSiblingElement("fracture");
	}
	stratum.sortFractures();
}