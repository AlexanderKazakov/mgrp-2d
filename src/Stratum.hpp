#ifndef STRATUM_HPP
#define	STRATUM_HPP

#include <vector>
#include <algorithm>
#include "Fracture.hpp"
#include "Field.hpp"
#include "Visualization.hpp"

class Stratum {
public:
	Stratum();
	~Stratum();
	void addFracture(const Fracture &fracture);
	void setRheology(double _G, double _nu);
	void sortFractures();
	int calculateNextFracture();
	Field calculateImpactInPoint(const double &x, const double &y);
	void visualize();
private:
	double G;	//	Shear modulus
	double nu;	//	Poisson's ratio
	double E;	//	Young's modulus	E = 2 * G * (1 + nu)
	std::vector<Fracture> fractures;
	std::vector<Fracture>::iterator currentFracture;
};

#endif	/* STRATUM_HPP */

