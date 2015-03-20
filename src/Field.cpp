#include "Field.hpp"

Field::Field() {
	Sxx = Sxy = Syy = 0;
// Uncomment this code if calculation of displacements become necessary
/*	Ux = Uy = 0; */
}

Field::~Field() {
}

void Field::operator+=(const Field &other) {
	Sxx += other.Sxx;
	Sxy += other.Sxy;
	Syy += other.Syy;
// Uncomment this code if calculation of displacements become necessary
/*	Ux  += other.Ux;
	Uy  += other.Uy; */
}

Field Field::inRotatedAxis(double beta) const {
	Field rotated;
// Uncomment this code if calculation of displacements become necessary
/*	rotated.Ux  = Ux * cos (beta) + Uy * sin (beta);
	rotated.Uy  = Uy * cos (beta) - Ux * sin (beta); */
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