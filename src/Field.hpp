#ifndef FIELD_HPP
#define	FIELD_HPP

#include <cmath>
#include <iostream>

class Field {
public:
	Field();
	~Field();
	void clear();
	void operator+=(const Field &other);
	Field inRotatedAxis(double beta);
	double directionOfMaxTensileStress();
	double Mises();
	double Sxx;
	double Syy;
	double Sxy;
	double Ux;
	double Uy;
private:

};

namespace std {
	inline std::ostream& operator<<(std::ostream &os, const Field &field) {
		os << field.Sxx << "\t" << field.Sxy << "\t" << field.Syy << std::endl;
		return os;
	};
}

#endif	/* FIELD_HPP */

