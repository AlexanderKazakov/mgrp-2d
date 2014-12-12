#ifndef STRATUM_HPP
#define	STRATUM_HPP

#include <vector>
#include "Fracture.hpp"
#include "Field.hpp"
#include <mgl2/mgl.h>

class Fracture;

class Stratum {
public:
	Stratum();
	~Stratum();
	void addFracture(int number, double x, double y, double beta, double h_length,
		int numOfBreaks, double a, double b, double c, std::string pressureType,
		std::string tip, std::string rotation);
	void setRheology(double _G, double _nu);
	void getRheology(double &_G, double &_nu);
	void setStresses(double _Sxx, double _Sxy, double _Syy);
	void setRanges(double _Xmin, double _Xmax, double _Ymin, double _Ymax);
	void reserve(int numberOfFractures);
	/**
	 * Calculate next one fracture supposing that all previous 
	 * fractures stay constant
     */
	void calculate();
	/**
	 * Calculate field in point (x, y) caused by existing fractures and
	 * external stresses
     * @param x	x-coord of point
     * @param y y-coord of point
     * @return Actual Field in (x, y)
     */
	Field calculateImpactInPoint(const double &x, const double &y);
	void visualize();
private:
	double G;	//	Shear modulus
	double nu;	//	Poisson's ratio
	double E;	//	Young's modulus	E = 2 * G * (1 + nu)
	double Xmin, Xmax, Ymin, Ymax;	//	Ranges for plot
	double Sxx;	//	
	double Sxy;	//	External stresses
	double Syy;	//
	std::vector<Fracture> fractures;
	std::vector<Fracture>::iterator currentFracture;
	void drawField(mglGraph &gr);
	void drawFractures(mglGraph &gr);
	void drawStressDirections(mglGraph &gr);
};

#endif	/* STRATUM_HPP */

