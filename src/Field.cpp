#include "Field.hpp"

Field::Field() {
	Sxx = Sxy = Syy = Ux = Uy = 0;
}

Field::~Field() {
}

void Field::operator+=(const Field &other) {
	Sxx += other.Sxx;
	Sxy += other.Sxy;
	Syy += other.Syy;
	Ux  += other.Ux;
	Uy  += other.Uy;
}

Field Field::inRotatedAxis(double beta) const {
	Field rotated;
	rotated.Ux  = Ux * cos (beta) + Uy * sin (beta);
	rotated.Uy  = Uy * cos (beta) - Ux * sin (beta);
	rotated.Sxx = Sxx * cos(beta) * cos(beta) 
	            + Sxy * sin (2 * beta) 
	            + Syy * sin(beta) * sin(beta);
	rotated.Syy = Syy * cos(beta) * cos(beta) 
	            - Sxy * sin (2 * beta) 
	            + Sxx * sin(beta) * sin(beta);
	rotated.Sxy = (Syy - Sxx) * sin(beta) * cos(beta) 
	            + Sxy * cos (2 * beta);
	return rotated;
}

double Field::directionOfMaxTensileStress() const {
	return arctan((Smax() - Sxx), Sxy);
}

double Field::Trace() const {
	return Sxx + Syy;
}

double Field::Smax() const {
	return (Sxx + Syy + sqrt((Sxx - Syy) * (Sxx - Syy) + 4 * Sxy * Sxy)) / 2;
}

double Field::arctan(const double& a, const double& b) const {
	double beta = atan(a / b);
	if (b == 0) {
		if (a > 0) beta = M_PI / 2;
		if (a < 0) beta = - M_PI / 2;
		if (a == 0) print("Error in atan: 0/0");
	}
	return beta;
}