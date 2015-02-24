#include "Breaker.hpp"

Breaker::Breaker() {
}

Breaker::~Breaker() {
}

void Breaker::setType(double _a, double _b, double _c, 
                      std::string _pressureType) {
	a = _a; b = _b; c = _c;
	pressureType = _pressureType;
}

void Breaker::getType(double& _a, double& _b, double& _c, 
                      std::string& _pressureType) const {
	_a = a; _b = b; _c = c;
	_pressureType = pressureType;
}

void Breaker::calculatePressure(std::vector<Element> &elementsL, 
                                std::vector<Element> &elementsR) const {
	
	if ( pressureType == "const" ) {
		for (auto element = elementsL.begin(); element != elementsL.end(); element++) {
			element->setSigmaN(c);
		}
		for (auto element = elementsR.begin(); element != elementsR.end(); element++) {
			element->setSigmaN(c);
		}
	} else if ( pressureType == "lag" ) {
		int i = 1;
		for (auto element = elementsL.begin(); element != elementsL.end(); element++) {
			double t = (i * 1.0) / elementsL.size();
			if (fabs(t) > a)
				element->setSigmaN(0);
			else
				element->setSigmaN(c);
			i++;
		}
		i = 1;
		for (auto element = elementsR.begin(); element != elementsR.end(); element++) {
			double t = (i * 1.0) / elementsR.size();
			if (fabs(t) > a)
				element->setSigmaN(0);
			else
				element->setSigmaN(c);
			i++;
		}
	} else if ( pressureType == "polynomial" ) {
		int i = 1;
		for (auto element = elementsL.begin(); element != elementsL.end(); element++) {
			double t = (i * 1.0) / elementsL.size();
			element->setSigmaN(a*t*t + b*t + c);
			i++;
		}
		i = 1;
		for (auto element = elementsR.begin(); element != elementsR.end(); element++) {
			double t = (i * 1.0) / elementsR.size();
			element->setSigmaN(a*t*t + b*t + c);
			i++;
		}
	}
	
}