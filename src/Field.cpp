#include "Field.hpp"

Field::Field() {
	Sxx = Sxy = Syy = Ux = Uy = 0;
}

Field::~Field() {
}

void Field::clear() {
	Sxx = Sxy = Syy = Ux = Uy = 0;
}

void Field::operator+=(const Field &other) {
	Sxx += other.Sxx;
	Sxy += other.Sxy;
	Syy += other.Syy;
	Ux += other.Ux;
	Uy += other.Uy;
}

Field Field::inRotatedAxis(double beta) {
	Field rotated;
	rotated.Ux = Ux * cos (beta) + Uy * sin (beta);
	rotated.Uy = Uy * cos (beta) - Ux * sin (beta);
	rotated.Sxx = Sxx * cos(beta) * cos(beta) + Sxy * sin (2 * beta) + Syy * sin(beta) * sin(beta);
	rotated.Syy = Syy * cos(beta) * cos(beta) - Sxy * sin (2 * beta) + Sxx * sin(beta) * sin(beta);
	rotated.Sxy = (Syy - Sxx) * sin(beta) * cos(beta) + Sxy * cos (2 * beta);
	return rotated;
}

double Field::directionOfMaxTensileStress() {
	return (atan((Smax() - Sxx) / Sxy));
}

double Field::Mises() {
	// TODO - Mises
	return 0;
}

double Field::Smax() {
	return (Sxx + Syy + sqrt((Sxx - Syy) * (Sxx - Syy) + 4 * Sxy * Sxy)) / 2;
}