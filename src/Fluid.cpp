#include "Fluid.hpp"

Fluid::Fluid() {
}

Fluid::~Fluid() {
}

void Fluid::calculatePressure(Break* breaks, int NumOfCalcBrks) {
	if ( (pressureType == "lag") || (pressureType == "const") ) {
		for (int i = -NumOfCalcBrks/2; i <= NumOfCalcBrks/2; i++) {
			breaks[i].setSigmaN(c);
		}
		if ( (pressureType == "lag") && (NumOfCalcBrks > 2) ) {
			breaks[-NumOfCalcBrks/2].setSigmaN(0);
			breaks[NumOfCalcBrks/2].setSigmaN(0);
		}
	}
	if ( pressureType == "polynomial" ) {
		for (int i = -NumOfCalcBrks/2; i <= NumOfCalcBrks/2; i++) {
			double t = (2.0 * i) / (NumOfCalcBrks - 1);
			breaks[i].setSigmaN(a*t*t + b*t + c);
		}
		if (NumOfCalcBrks == 1)
			breaks[0].setSigmaN(c);
	}

}

void Fluid::setType(double _a, double _b, double _c, std::string _pressureType) {
	a = _a;	b = _b; c = _c;
	pressureType = _pressureType;
}