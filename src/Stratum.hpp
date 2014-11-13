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
	void setRanges(double _Xmin, double _Xmax, double _Ymin, double _Ymax);
	void sortFractures();
	int calculateNextFracture();
	Field calculateImpactInPoint(const double &x, const double &y);
	void visualize();
private:
	double G;	//	Shear modulus
	double nu;	//	Poisson's ratio
	double E;	//	Young's modulus	E = 2 * G * (1 + nu)
	double Xmin, Xmax, Ymin, Ymax;	//	Ranges for plot
	std::vector<Fracture> fractures;
	std::vector<Fracture>::iterator currentFracture;
	void drawField(mglGraph &gr);
	void drawFractures(mglGraph &gr);
	void drawDirections(mglGraph &gr);
};

#endif	/* STRATUM_HPP */

