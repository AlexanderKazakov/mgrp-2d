#ifndef FLUID_HPP
#define	FLUID_HPP

#include <string>
#include <vector>
#include "Break.hpp"

class Fluid {
public:
	Fluid();
	~Fluid();
	void setType(double _a, double _b, double _c, std::string _pressureType);
	void setPressure(std::vector<Break> &breaks);
private:
	double a, b, c;
	std::string pressureType;
};

#endif	/* FLUID_HPP */

