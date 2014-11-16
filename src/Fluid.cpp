#include "Fluid.hpp"

Fluid::Fluid() {
}

Fluid::~Fluid() {
}

void Fluid::setPressure(std::vector<Break> &breaks) {
	std::vector<Break>::iterator break1 = breaks.begin();	
	if ( (pressureType == "lag") || (pressureType == "const") ) {
		while (break1 != breaks.end()) {
			break1->setSigmaN(c);
			break1++;
		}
		if ( (pressureType == "lag") && (breaks.size() > 2) ) {
			breaks.front().setSigmaN(0);
			breaks.back().setSigmaN(0);
		}
	}
	if ( pressureType == "polynomial" ) {
		while (break1 != breaks.end()) {
			double t = (2.0 * break1->getNumber()) / (breaks.size() - 1);
			break1->setSigmaN(a*t*t + b*t + c);
			//std::cout << a*t*t + b*t + c << "\t";
			break1++;
		}
		//std::cout << std::endl;
		if (breaks.size() == 1)
			breaks.front().setSigmaN(c);
	}
}

void Fluid::setType(double _a, double _b, double _c, std::string _pressureType) {
	a = _a;	b = _b; c = _c;
	pressureType = _pressureType;
}