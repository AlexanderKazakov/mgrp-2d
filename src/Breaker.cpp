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

void Breaker::calculatePressure(Element* lmnts, int numOfCalcLmnts) const {
	if ( pressureType == "const" ) {
		for (int i = -numOfCalcLmnts/2; i <= numOfCalcLmnts/2; i++)
			lmnts[i].setSigmaN(c);
	}
	
	if ( pressureType == "lag" ) {
		for (int i = -numOfCalcLmnts/2; i <= numOfCalcLmnts/2; i++) {
			double t = (2.0 * i) / (numOfCalcLmnts - 1);
			if ( fabs(t) > a )
				lmnts[i].setSigmaN(0);
			else
				lmnts[i].setSigmaN(c);
		}
	}
	
	if ( pressureType == "polynomial" ) {
		for (int i = -numOfCalcLmnts/2; i <= numOfCalcLmnts/2; i++) {
			double t = (2.0 * i) / (numOfCalcLmnts - 1);
			lmnts[i].setSigmaN(a*t*t + b*t + c);
		}
		if (numOfCalcLmnts == 1)
			lmnts[0].setSigmaN(c);
	}

}