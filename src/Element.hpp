#ifndef ELEMENT_HPP
#define	ELEMENT_HPP

#include <cmath>
#include <iostream>

#include "Field.hpp"
#include "util.hpp"

/**
 * Class that represent the boundary element in Crouch's 
 * displacement discontinuity method for calculation fractures.
 * 
 */

class Element {
public:
	double Ds; // displacement discontinuity along element
	double Dn; // displacement discontinuity across element
	double Cx; // Center on x
	double Cy; // Center on y
	double beta; //	Angle to x
	
	Element();
	/**
	 * Constructor for boundary element
     * @param a	half-length	
     * @param Cx x-coord of the center of the element
     * @param Cy y-coord of the center of the element
     * @param beta angle to x-axis
     * @param G shear modulus of stratum
     * @param nu Poisson's ratio of stratum
     */
	Element(double a, double Cx, double Cy, double beta, double G, double nu);
	~Element();
	/**
	 * Get the right-hand side corresponding to Ds of this element 
	 * in the linear system on displacement discontinuities 
     * @return right-hand side of SLE corresponding to Ds
     */
	double getBs() const;
	/**
	 * Get the right-hand side corresponding to Dn of this element 
	 * in the linear system on displacement discontinuities 
     * @return right-hand side of SLE corresponding to Dn
     */
	double getBn() const;
	/**
     * @return half-length of the element
     */
	double getA() const;
	/**
     * @param _sigmaN value to set as sigmaN
     */
	void setSigmaN(const double &_sigmaN);
	/**
	 * Set impact of external fractures and stresses at infinity on the element
     * @param field 
     */
	void setExternalImpact(const Field &field);
	/**
	 * Calculate the components of the matrix of linear system on displacement 
	 * discontinuities corresponding to impact of this element on element lmnt2
     * @param lmnt2 element to calculate impact on
     * @param Ass impact of this->Ds on lmnt2.Ds
     * @param Asn impact of this->Ds on lmnt2.Dn
     * @param Ans impact of this->Dn on lmnt2.Ds
     * @param Ann impact of this->Dn on lmnt2.Dn
     */
	void calculateImpactOn(const Element &lmnt2, double &Ass, double &Asn,
	                                             double &Ans, double &Ann) const;
	/**
	 * Calculate impact of this element on field in the point (x_glob, y_glob)
     * @param x_glob x-coord in global system of coordinates
     * @param y_glob y-coord in global system of coordinates
     * @return field in (x_glob, y_glob)
     */
	Field calculateImpactInPoint(const double &x_glob, 
	                             const double &y_glob) const;
	/**
     * @return Mode I stress intensity factor near the element
     */
	double K1() const;
	/**
     * @return Mode II stress intensity factor near the element
     */
	double K2() const;
	
private:
	double G, nu; // Rheology parameters
	double a; // Half-length
	double sigmaN; // Impact of the inner breaker
	double sigmaS; // Impact of the inner breaker
	double externalSigmaN; // Impact of already existing fractures 
	double externalSigmaS; // Impact of already existing fractures
	
	// Functions used in calculation of impact of element in point (x, y)
	double F1(const double &x, const double &y) const;
	double F2(const double &x, const double &y) const;
	double F3(const double &x, const double &y) const;
	double F4(const double &x, const double &y) const;
	double F5(const double &x, const double &y) const;
	double F6(const double &x, const double &y) const;
	double F7(const double &x, const double &y) const;
	
};

namespace std {
	inline std::ostream& operator<<(std::ostream &os, const Element &brk) {
		os << "\t" << brk.Cx
				//<< "\t" << brk.Ds << "\t" << brk.Dn 
				<< "\t" << brk.getBs() << "\t" << brk.getBn() 
				<< std::endl;
		return os;
	};
}

#endif	/* ELEMENT_HPP */

