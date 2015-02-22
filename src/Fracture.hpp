#ifndef FRACTURE_HPP
#define	FRACTURE_HPP

#include <string>
#include <iostream>

#include <gsl/gsl_linalg.h>

#include "Element.hpp"
#include "Field.hpp"
#include "Breaker.hpp"
#include "Stratum.hpp"
#include "util.hpp"


class Stratum;

/**
 * Class that contains displacement discontinuity boundary elements (lmnts)
 * and acts on them when fracture is growing.
 *
 */

class Fracture {
public:
	Fracture();
	/**
	 * Constructor of the fracture.
	 * Pressure of the inner breaker of the fracture is equal to
	 * - (at^2 + bt + c), t from -1 to 1, if pressureType is polynomial,
	 * and ( -c ), if pressureType is const
     * @param stratum pointer to the stratum which owns the fracture
     * @param number index number of fracture
     * @param halfLengthOfLmnts half-length of boundary elements of fracture
     * @param numOfLmntsL required number of boundary elements at the 
	 * left part at the end of computation
	 * @param numOfLmntsR required number of boundary elements at the 
	 * right part at the end of computation
     * @param _a pressure parameter
     * @param _b pressure parameter
     * @param _c pressure parameter 
     * @param pressureType type of pressure of the inner breaker
     * @param tip rule to use special boundary element at the tip of the
	 * fracture, "const" or TODO
     * @param rotation rule to use splitting on actual configuration's
	 * calculation and calculation of direction of the fracture's growth,
	 * "predictor" or "predictor-corrector"
     */
	Fracture(Stratum *stratum, int number, int numOfLmntsL, int numOfLmntsR, 
	         double halfLengthOfLmnts, double _a, double _b, double _c,
	         std::string pressureType, std::string tip, std::string rotation);
	~Fracture();
	/**
	 * Dynamically (!) allocate memory for lmnts and create the initial lmnt
     * @param x x-coord of the initial lmnt
     * @param y y-coord of the initial lmnt
     * @param beta angle to x of the initial lmnt
     */
	void allocateLmnts(double x, double y, double beta);
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
	 * Get number of elements on the left edge of the fracture
     * @return number of elements on the left edge of the fracture
     */
	int getNumOfLmntsL() const;
	/**
	 * Get number of elements on the right edge of the fracture
     * @return number of elements on the right edge of the fracture
     */
	int getNumOfLmntsR() const;
	/**
	 * Get points for drawing the fracture on the graph
     * @param x
     * @param y
     */
	void getPointsForPlot(double *x, double *y) const;
	/**
	 * Get displacements v(x) for drawing on the graph 
	 * (for testing purposes)
     * @param x
     * @param v
     */
	void getPointsForDisplacementPlot(double* x, double* v) const;
	/**
     * @return length of the left part of the fracture
     */
	double getHalfLengthL() const;
	/**
     * @return length of the right part of the fracture
     */
	double getHalfLengthR() const;
	
	// It's necessary for fast and easy setting the pressure of the inner 
	// breaker on every step of calculation.
	friend class Breaker;
private:
	// if the fraction is stopped by compression on the left edge
	bool fractionIsStoppedL;
	// if the fraction is stopped by compression on the right edge
	bool fractionIsStoppedR;
	int number; // index number of fracture
	// number of required elements on the left edge
	int numOfLmntsL; 
	// number of required elements on the right edge
	int numOfLmntsR; 
	double G, nu; // rheology parameters of stratum
	double a; // half-length of boundary elements
	double halfLengthL; // length of the left part of the fracture
	double halfLengthR; // length of the right part of the fracture
	std::string tip; // rule to use special boundary element at the tip
	// rule to use splitting on actual configuration's calculation 
	// and calculation of direction of the fracture's growth
	std::string rotation;
		
	// displacement discontinuity boundary elements 
	// at the left part of the fracture
	std::vector<Element> lmntsL;	
	// displacement discontinuity boundary elements 
	// at the right part of the fracture
	std::vector<Element> lmntsR;	
	Breaker breaker;	//	class of the inner breaker
	Stratum *stratum;	//	pointer to the stratum which owns the fracture
	/**
	 * Fill in the system of linear equations on values of displacement 
	 * discontinuity of elements, solve it by GSL library and set new values
	 * of displacement discontinuity to elements.
	 * 
	 * It's the main idea of the program. It uses displacement discontinuity 
	 * boundary method proposed by Crouch.
     */
	void calculateLmnts();
	/**
	 * Add new elements to the fracture. Three small tip elements in both edges 
	 * are replaced by one ordinary element with equal length and angle. Next, 
	 * new three small elements are added to the edges with new angle from parameters.
     * @param deltaBeta1 rotation of the left tip element towards its neighbour
     * @param deltaBeta2 rotation of the right tip element towards its neighbour
     */
	void grow(const double &deltaBeta1, const double &deltaBeta2);
	/**
	 * Add two elements to the edges of fracture.
     * @param deltaBeta1 rotation of the first (left) element towards its neighbour
     * @param deltaBeta2 rotation of the last (right) element towards its neighbour
	 * @param half_length half-length of the elements to add
     */
	void addNewLmnts(const double &deltaBeta1, const double &deltaBeta2,
	                                            const double &half_length);
	/**
	 * Replace three small tip elements on both edges with new three small tip 
	 * elements with other angle from parameters
     * @param deltaBeta1 rotation of the left tip element towards its neighbour
     * @param deltaBeta2 rotation of the right tip element towards its neighbour
     */
	void replaceTipElements(const double &deltaBeta1, const double &deltaBeta2);
	/**
	 * Calculate the angle of rotation of the fracture propagation
	 * near the given lmnt 
     * @param lmnt1 lmnt to calculate the rotation of propagation direction
     * @return angle of rotation of propagation direction
     */
	double calcAngleOfRotation(const Element &lmnt1) const;
};

#endif	/* FRACTURE_HPP */

