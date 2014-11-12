#ifndef STRATUM_HPP
#define	STRATUM_HPP

#include <vector>
#include <algorithm>
#include "Fracture.hpp"
#include "Field.hpp"
#include <mgl2/mgl.h>

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
	void drawField(mglGraph &gr, const double &Xmin, const double &Xmax,
									const double &Ymin, const double &Ymax);
	void drawFractures(mglGraph &gr);
	void drawDirections(mglGraph &gr, const double &Xmin, const double &Xmax,
									const double &Ymin, const double &Ymax);
};

#endif	/* STRATUM_HPP */

