#ifndef FRACTURE_HPP
#define	FRACTURE_HPP

#include <string>
#include <iostream>
#include <gsl/gsl_linalg.h>
#include "Break.hpp"
#include "Field.hpp"
#include "Fluid.hpp"
#include "Stratum.hpp"

class Stratum;

class Fracture {
public:
	Fracture();
	Fracture(Stratum *stratum, int number, double h_length,
		int numOfBreaks, double a, double b, double c, std::string pressureType);
	~Fracture();
	void allocateBreaks(double x, double y, double beta);
	int getNumber() const;
	bool operator==(const Fracture &other) const;
	bool operator!=(const Fracture &other) const;
	bool operator<(const Fracture &other) const;
	void calculate();
	Field calculateImpactInPoint (const double &x, const double &y) const;
	int getNumOfBreaks() const;
	void getPointsForPlot(double *x, double *y) const;
	friend class Fluid;
private:
	int number;
	int middle;
	int front;
	int back;
	int numOfBrks; //	number of elements in the fracture
	int numOfCalcBrks;	//	number of already calculated elements
	double G, nu;
	double half_lengthOfBreaks;
	Break *breaks;
	Fluid fluid;
	Stratum *stratum;
	bool calculateBreaks();
	double calcAngleOfRotation(const double &K1, const double &K2) const;
};

#endif	/* FRACTURE_HPP */

