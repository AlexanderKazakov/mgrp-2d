#ifndef FIELD_HPP
#define	FIELD_HPP

#include <cmath>
#include <iostream>

#include "util.hpp"


/**
 * Class for some auxiliary actions with linear elastic field 
 * in some point (x, y)
 */


class Field {
public:
	Field();
	~Field();
	/**
	 * Add other field to this field
     * @param other field to add
     */
	void operator+=(const Field &other);
	/**
	 * Calculate field in the same point in rotated axis
     * @param beta rotation angle of new axis to old axis
     * @return field in rotated axis
     */
	Field inRotatedAxis(double beta) const;
	/**
     * @return angle to x of direction of maximal tensile stresses
     */
	double directionOfMaxTensileStress() const;
	/**
     * @return Sxx + Syy 
     */
	double Trace() const;
	/**
     * @return maximal eigenvalue of stress tensor 
     */
	double Smax() const;
	double Sxx;
	double Syy;
	double Sxy;
// Uncomment this code if calculation of displacements become necessary
/*	double Ux;
	double Uy; */
private:

};

namespace std {
	inline std::ostream& operator<<(std::ostream &os, const Field &field) {
		os << field.Sxx << "\t" << field.Sxy << "\t" << field.Syy << std::endl;
		return os;
	};
}

#endif	/* FIELD_HPP */

