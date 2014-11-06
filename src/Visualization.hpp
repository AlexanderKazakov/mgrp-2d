#ifndef VISUALIZATION_HPP
#define	VISUALIZATION_HPP

#include <mgl2/mgl.h>
#include <iostream>
#include "Field.hpp"

class Visualization {
public:
	Visualization();
	~Visualization();
	void sample1();
	void sample2();
	int densPlot(mglGraph *gr);
	void plotFracture(int N, double *x, double *y);
	//void plotField(Field Stratum::(*f)(const double&, const double&));

private:
	mglGraph gr;
};

#endif	/* VISUALIZATION_HPP */

