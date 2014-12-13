#ifndef FRACTURE_HPP
#define	FRACTURE_HPP

#include <string>
#include <iostream>
#include <gsl/gsl_linalg.h>
#include "Break.hpp"
#include "Field.hpp"
#include "Fluid.hpp"
#include "Stratum.hpp"
#include "util.hpp"


class Stratum;

/**
 * Class that contains displacement discontinuity boundary elements (breaks)
 * and acts on them when fracture is growing.
 *
 */

class Fracture {
public:
	Fracture();
	/**
	 * Constructor of the fracture.
	 * Pressure of the inner fluid of the fracture is equal to
	 * - (at^2 + bt + c), t from -1 to 1, if pressureType is polynomial,
	 * and ( -c ), if pressureType is const
     * @param stratum pointer to the stratum which owns the fracture
     * @param number index number of fracture
     * @param h_length half-length of boundary elements of fracture
     * @param numOfBreaks required number of boundary elements at the end of
	 * computation
     * @param a pressure parameter
     * @param b pressure parameter
     * @param c pressure parameter 
     * @param pressureType type of pressure of the inner fluid
     * @param tip rule to use special boundary element at the tip of the
	 * fracture, "const" or "sqrt" 
     * @param rotation rule to use splitting on actual configuration's
	 * calculation and calculation of direction of the fracture's growth,
	 * "predictor" or "predictor-corrector"
     */
	Fracture(Stratum *stratum, int number, double h_length,	int numOfBreaks,
			double a, double b, double c, std::string pressureType, 
			std::string tip, std::string rotation);
	~Fracture();
	/**
	 * Dynamically (!) allocate memory for breaks and create the initial break
     * @param x x-coord of the initial break
     * @param y y-coord of the initial break
     * @param beta angle to x of the initial break
     */
	void allocateBreaks(double x, double y, double beta);
	/**
	 * Calculate itself taking into account already existing in the stratum 
	 * stresses. Do calculation until the number of elements is less than 
	 * required and the fracture is not stopped by compression of its tips 
     */
	void calculate();
	/**
	 * Calculate field in point (x, y) caused by the fracture
     * @param x
     * @param y
     * @return field in (x, y) caused by the fracture
     */
	Field calculateImpactInPoint (const double &x, const double &y) const;
	/**
	 * Get number of elements in the fracture
     * @return number of elements in the fracture
     */
	int getNumOfBreaks() const;
	/**
	 * Get points for drawing the fracture on the graph
     * @param x
     * @param y
     */
	void getPointsForPlot(double *x, double *y) const;
	// It's necessary for fast and easy setting the pressure of the inner 
	// fluid on every step of calculation.
	friend class Fluid;
private:
	bool fractionIsStopped;	//	if the fraction is stopped 
							//	by compression on the edges
	int number;	//	index number of fracture
	int middle;	//	index of the initial element of the fracture
	int front;	//	index of the one tip element of the fracture
	int back;	//	index of the another tip element of the fracture
	int numOfBrks; //	number of required elements
	int numOfCalcBrks;	//	number of already calculated elements
	double G, nu;	//	rheology parameters
	double half_lengthOfBreaks;	//	half-length of boundary elements
	std::string tip; //	rule to use special boundary element at the tip
	std::string rotation;
		//	rule to use splitting on actual configuration's calculation 
		//	and calculation of direction of the fracture's growth
	Break *breaks;	//	displacement discontinuity boundary elements
	Fluid fluid;	//	class of the inner fluid
	Stratum *stratum;	//	pointer to the stratum which owns the fracture
	/**
	 * Fill in the system of linear equations on values of displacement 
	 * discontinuity of elements, solve it by GSL library and set new values
	 * of displacement discontinuity to elements.
	 * 
	 * It's the main idea of the program. It uses displacement discontinuity 
	 * boundary method proposed by Crouch.
     */
	void calculateBreaks();
	/**
	 * Add two tip elements to the edges of fracture
     * @param deltaBeta1 rotation of the first element towards its neighbour
     * @param deltaBeta2 rotation of the last element towards its neighbour
     */
	void addNewBreaks(const double &deltaBeta1, const double &deltaBeta2);
	/**
	 * Calculate the angle of rotation of the fracture propagation
	 * near the given break 
     * @param break1 break to calculate the rotation of propagation direction
     * @return angle of rotation of propagation direction
     */
	double calcAngleOfRotation(const Break &break1) const;
};

#endif	/* FRACTURE_HPP */

