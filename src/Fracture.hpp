#ifndef FRACTURE_HPP
#define	FRACTURE_HPP

#include <vector>
#include <string>
#include <iostream>
#include <gsl/gsl_linalg.h>
#include "Break.hpp"
#include "Field.hpp"
#include "Fluid.hpp"

class Fracture {
public:
	Fracture();
	Fracture(int number, double x, double y, double beta, double length,
		int numOfBreaks, double a, double b, double c, double G, double nu,
													std::string pressureType);
	~Fracture();
	int getNumber() const;
	bool operator==(const Fracture &other) const;
	bool operator!=(const Fracture &other) const;
	bool operator<(const Fracture &other) const;
	void calculate(std::vector<Fracture>::const_iterator firstFracture);
	Field calculateImpactInPoint (const double &x, const double &y) const;
	int getNumOfBreaks() const;
	void getPointsForPlot(double *x, double *y) const;
	friend class Fluid;
private:
	int number;
	int numOfBreaks;
	int numOfCalculatedBreaks;
	double G, nu;
	double half_lengthOfBreaks;
	std::vector<Break> breaks;
	Fluid fluid;
	bool calculateBreaks();
	double calcAngleOfRotation(const double &K1, const double &K2) const;
};

#endif	/* FRACTURE_HPP */

