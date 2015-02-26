#include "Breaker.hpp"
#include "Fracture.hpp"

Breaker::Breaker() {
	leftN = rightN = 0;
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

void Breaker::calculatePressure(Fracture *frac) {
	if (frac->breakerIsInjected) {
		
		if ( pressureType == "const" ) {
			for (auto element = frac->elementsL.begin(); 
			          element != frac->elementsL.end(); element++) {
				element->setSigmaN(c);
			}
			for (auto element = frac->elementsR.begin(); 
			          element != frac->elementsR.end(); element++) {
				element->setSigmaN(c);
			}
		} else if ( pressureType == "lag" ) {
			int i = 1;
			for (auto element = frac->elementsL.begin(); 
			          element != frac->elementsL.end(); element++) {
				double t = (i * 1.0) / frac->elementsL.size();
				if (fabs(t) > a)
					element->setSigmaN(0);
				else
					element->setSigmaN(c);
				i++;
			}
			i = 1;
			for (auto element = frac->elementsR.begin(); 
			          element != frac->elementsR.end(); element++) {
				double t = (i * 1.0) / frac->elementsR.size();
				if (fabs(t) > a)
					element->setSigmaN(0);
				else
					element->setSigmaN(c);
				i++;
			}
		} else if ( pressureType == "polynomial" ) {
			int i = 1;
			for (auto element = frac->elementsL.begin(); 
			          element != frac->elementsL.end(); element++) {
				double t = (i * 1.0) / frac->elementsL.size();
				element->setSigmaN(a*t*t + b*t + c);
				i++;
			}
			i = 1;
			for (auto element = frac->elementsR.begin(); 
			          element != frac->elementsR.end(); element++) {
				double t = (i * 1.0) / frac->elementsR.size();
				element->setSigmaN(a*t*t + b*t + c);
				i++;
			}
		}
		
	} else {
		if ( ! leftN ) {
			leftN = frac->elementsL.size();
			rightN = frac->elementsR.size();
		}
		
		if ( pressureType == "const" ) {
			int i = 1;
			for (auto element = frac->elementsL.begin(); 
			          element != frac->elementsL.end(); element++) {
				if (i < leftN) element->setSigmaN(eim * c);
				else element->setSigmaN(0);
				i++;
			}
			i = 1;
			for (auto element = frac->elementsR.begin(); 
			          element != frac->elementsR.end(); element++) {
				if (i < rightN) element->setSigmaN(eim * c);
				else element->setSigmaN(0);
				i++;
			}
		} else if ( pressureType == "lag" ) {
			int i = 1;
			for (auto element = frac->elementsL.begin(); 
			          element != frac->elementsL.end(); element++) {
				double t = (i * 1.0) / leftN;
				if (fabs(t) > a)
					element->setSigmaN(0);
				else
					element->setSigmaN(eim * c);
				i++;
			}
			i = 1;
			for (auto element = frac->elementsR.begin(); 
			          element != frac->elementsR.end(); element++) {
				double t = (i * 1.0) / rightN;
				if (fabs(t) > a)
					element->setSigmaN(0);
				else
					element->setSigmaN(eim * c);
				i++;
			}
		} else if ( pressureType == "polynomial" ) {
			int i = 1;
			for (auto element = frac->elementsL.begin(); 
			          element != frac->elementsL.end(); element++) {
				double t = (i * 1.0) / leftN;
				if (fabs(t) > 1) element->setSigmaN(0);
				else element->setSigmaN(eim * (a*t*t + b*t + c));
				i++;
			}
			i = 1;
			for (auto element = frac->elementsR.begin(); 
			          element != frac->elementsR.end(); element++) {
				double t = (i * 1.0) / rightN;
				if (fabs(t) > 1) element->setSigmaN(0);
				else element->setSigmaN(eim * (a*t*t + b*t + c));
				i++;
			}
		}
		
	}
}