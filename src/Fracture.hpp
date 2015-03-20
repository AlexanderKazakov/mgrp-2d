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
class Breaker;

/**
 * Class that contains displacement discontinuity boundary elements (elements)
 * and acts on them when fracture is growing.
 *
 */

class Fracture {
private:
	// if the fraction is stopped by compression on the left corner
	bool fractionIsStoppedL;
	// if the fraction is stopped by compression on the right corner
	bool fractionIsStoppedR;
	int number; // index number of fracture
	double a; // half-length of boundary elements
	// actual area of the fracture (smth like length * displacementDiscontinuity)
	double volume;
	// area of the fracture (smth like length * displacementDiscontinuity),
	// required in task
	double taskVolume;
	double leftLength; // length of the left part of the fracture
	double rightLength; // length of the right part of the fracture
	// change of the direction of fracture's propagation
	// on predictor step of predictor-corrector method, it's necessary
	// on corrector step, so we keep it here
	double deltaBetaL, deltaBetaR;
	static std::string tip; // rule to use special boundary elements at the tip
	
	// displacement discontinuity boundary elements 
	// at the left part of the fracture
	std::vector<Element> elementsL;	
	// displacement discontinuity boundary elements 
	// at the right part of the fracture
	std::vector<Element> elementsR;	
	Breaker *breaker; // pointer to the class of the inner breaker
	static Stratum *stratum; // pointer to the stratum which owns the fracture
public:
	// if the breaker is injecting to the fracture now, used by breaker 
	// to set sigmaN to elements
	bool breakerIsInjecting;
	Fracture();
	/**
	 * Constructor of the fracture.
     * @param number index number of fracture
     * @param volume required in task volume of the fracture at the end of injection 
     * @param halfLengthOfElements half-length of boundary elements of fracture
     */
	Fracture(int number, double volume, double halfLengthOfElements);
	~Fracture();
	/**
     * @param x x-coord of the initial element
     * @param y y-coord of the initial element
     * @param beta angle to x of the initial element
	 * @param _breaker Breaker with parameters to set for breaker of this fracture 
     */
	void allocateElements(double x, double y, double beta, Breaker _breaker);
	/**
     * @return true if required in task volume has been reached
	 * or fracture is stopped on both corners
     */
	bool isCompleted() const;
	/**
     * @return true if fracture is stopped on both corners
     */
	bool isStopped() const;
	/**
	 * Calculate field in point (x, y) caused by the fracture
     * @param x
     * @param y
     * @return field in (x, y) caused by the fracture
     */
	Field calculateImpactInPoint(const double &x, const double &y) const;
	/**
	 * Add new elements to the fracture. Three small tip elements in both corners 
	 * are replaced by one ordinary element with equal length and angle. Next, 
	 * new three small elements are added to the corners with new angle 
     */
	void grow();
	/**
	 * Replace three small tip elements on both corners with new three small tip 
	 * elements with other angle from predictor-corrector method
     */
	void correctRotation();
	/**
	 * Set external impact and pressure of breaker for all
	 * elements in the fracture.
     */
	void setExternalImpactAndBreakerPressure();
	/**
	 * Fill in the corresponding rectangular part of the matrix of SLE on
	 * A * x = b
	 * displacement discontinuities.
	 * The part starts from i1'th column, i2'th string of A.  
     * @param frac2 fracture which impact on this fracture will be calculated
     * @param A matrix of SLE
     * @param i1 column to start from
     * @param i2 string to start from
     */
	void fillInMatrixA(const Fracture &frac2, gsl_matrix *A,
	                              int i1, int i2) const;
	/**
	 * Fill in the corresponding part of the right-hand-side vector of SLE
	 * A * x = b
	 * on displacement discontinuities.
     * @param b vector to fill in
     * @param i index to start from
     */
	void fillInVectorB(gsl_vector *b, int i) const;
	/**
	 * Take values of displacement discontinuities from corresponding part of 
	 * vector x of SLE
	 * A * x = b
     * @param x vector to take from
     * @param i index to start from
     */
	void takeDDfromVectorX(const gsl_vector *x, int i);
	/**
     * @return index number of the fracture
     */
	int getNumber() const;
	/**
	 * Get number of elements in the fracture
     * @return number of elements in the fracture
     */
	int getNumOfElements() const;
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
	double getLeftLength() const;
	/**
     * @return length of the right part of the fracture
     */
	double getRightLength() const;
	/**
     * @return pointer to its breaker
     */
	Breaker *getBreaker() const;
	/**
	 * Set static parameter - pointer to stratum that contains all fractures
     * @param _stratum
     */
	void setStratum(Stratum *_stratum);
	/**
	 * Set static parameters - rule to use special boundary element at the corners
     * @param _tip
     */
	void setTip(std::string _tip);
	
	// It's necessary for fast and easy setting the pressure of the inner 
	// breaker on every step of calculation.
	friend class Breaker;

private:
	/**
	 * Add new elements to the fracture. Three small tip elements in both corners 
	 * are replaced by one ordinary element with equal length and angle. Next, 
	 * new three small elements are added to the corners with new angle from parameters.
     * @param deltaBeta1 rotation of the left tip element towards its neighbour
     * @param deltaBeta2 rotation of the right tip element towards its neighbour
     */
	void grow(const double &deltaBeta1, const double &deltaBeta2);
	/**
	 * Add two elements to the corners of fracture.
     * @param deltaBeta1 rotation of the first (left) element towards its neighbour
     * @param deltaBeta2 rotation of the last (right) element towards its neighbour
	 * @param half_length half-length of the elements to add
     */
	void addNewElements(const double &deltaBeta1, const double &deltaBeta2,
	                                            const double &half_length);
	/**
	 * Replace three small tip elements on both corners with new three small tip 
	 * elements with other angle from parameters
     * @param deltaBeta1 rotation of the left tip element towards its neighbour
     * @param deltaBeta2 rotation of the right tip element towards its neighbour
     */
	void replaceTipElements(const double &deltaBeta1, const double &deltaBeta2);
};

namespace std {
	inline std::ostream& operator<<(std::ostream &os, const Fracture &frac) {
		os << frac.getNumber()
				<< std::endl;
		return os;
	};
}

#endif	/* FRACTURE_HPP */

