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
	Fracture(Stratum *stratum, int number, double h_length,	int numOfBreaks,
			double a, double b, double c, std::string pressureType, 
			std::string tip, std::string rotation);
	~Fracture();
	void allocateBreaks(double x, double y, double beta);
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
	std::string tip; 
	std::string rotation;
	Break *breaks;
	Fluid fluid;
	Stratum *stratum;
	bool calculateBreaks();
	void addNewBreaks(const double &deltaBeta1, const double &deltaBeta2);
	double calcAngleOfRotation(const Break &break1) const;
};

#endif	/* FRACTURE_HPP */

