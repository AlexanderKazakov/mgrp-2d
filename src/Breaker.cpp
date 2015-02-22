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

void Breaker::calculatePressure(std::vector<Element> &lmntsL, 
                                std::vector<Element> &lmntsR) const {
	
	if ( pressureType == "const" ) {
		for (auto lmnt = lmntsL.begin(); lmnt != lmntsL.end(); lmnt++) {
			lmnt->setSigmaN(c);
		}
		for (auto lmnt = lmntsR.begin(); lmnt != lmntsR.end(); lmnt++) {
			lmnt->setSigmaN(c);
		}
	} else if ( pressureType == "lag" ) {
		int i = 1;
		for (auto lmnt = lmntsL.begin(); lmnt != lmntsL.end(); lmnt++) {
			double t = (i * 1.0) / lmntsL.size();
			if (fabs(t) > a)
				lmnt->setSigmaN(0);
			else
				lmnt->setSigmaN(c);
			i++;
		}
		i = 1;
		for (auto lmnt = lmntsR.begin(); lmnt != lmntsR.end(); lmnt++) {
			double t = (i * 1.0) / lmntsR.size();
			if (fabs(t) > a)
				lmnt->setSigmaN(0);
			else
				lmnt->setSigmaN(c);
			i++;
		}
	} else if ( pressureType == "polynomial" ) {
		int i = 1;
		for (auto lmnt = lmntsL.begin(); lmnt != lmntsL.end(); lmnt++) {
			double t = (i * 1.0) / lmntsL.size();
			lmnt->setSigmaN(a*t*t + b*t + c);
			i++;
		}
		i = 1;
		for (auto lmnt = lmntsR.begin(); lmnt != lmntsR.end(); lmnt++) {
			double t = (i * 1.0) / lmntsR.size();
			lmnt->setSigmaN(a*t*t + b*t + c);
			i++;
		}
	}
	
}