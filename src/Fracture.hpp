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
 * Class that contains displacement discontinuity boundary elements (elements)
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
     * @param halfLengthOfElements half-length of boundary elements of fracture
     * @param numOfElements required number of boundary elements at the 
	 * end of computation
     * @param _a pressure parameter
     * @param _b pressure parameter
     * @param _c pressure parameter 
     * @param pressureType type of pressure of the inner breaker
     * @param tip rule to use special boundary element at the tip of the
	 * fracture, "const" or TODO
     */
	Fracture(Stratum *stratum, int number, int numOfElements,
	         double halfLengthOfElements, double _a, double _b, double _c,
	         std::string pressureType, std::string tip);
	~Fracture();
	/**
     * @param x x-coord of the initial element
     * @param y y-coord of the initial element
     * @param beta angle to x of the initial element
     */
	void allocateElements(double x, double y, double beta);
	/**
     * @return true if required in task number of elements has been calculated
	 * or fracture is stopped on both corners
     */
	bool isCompleted() const;
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
	 * Replace three small tip elements on both edges with new three small tip 
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
	 * Get number of elements at the left part of the fracture
     * @return number of elements at the left part of the fracture
     */
	int getNumOfElementsL() const;
	/**
	 * Get number of elements at the right part of the fracture
     * @return number of elements at the right part of the fracture
     */
	int getNumOfElementsR() const;
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
	
	// It's necessary for fast and easy setting the pressure of the inner 
	// breaker on every step of calculation.
	friend class Breaker;
private:
	// if the fraction is stopped by compression on the left edge
	bool fractionIsStoppedL;
	// if the fraction is stopped by compression on the right edge
	bool fractionIsStoppedR;
	int number; // index number of fracture
	// number of required elements at the end of computation
	int numOfElements;
	double G, nu; // rheology parameters of stratum
	double a; // half-length of boundary elements
	double leftLength; // length of the left part of the fracture
	double rightLength; // length of the right part of the fracture
	// change of the direction of fracture's propagation
	// on predictor step of predictor-corrector method, it's necessary
	// on corrector step, so we keep it here
	double deltaBetaL, deltaBetaR;
	std::string tip; // rule to use special boundary element at the tip
	
	// displacement discontinuity boundary elements 
	// at the left part of the fracture
	std::vector<Element> elementsL;	
	// displacement discontinuity boundary elements 
	// at the right part of the fracture
	std::vector<Element> elementsR;	
	
//	// class to iterate all the elements of the fracture by one cycle
//	class FracIter {
//	private:
//		//std::vector<Element> *elementsL;
//		std::vector<Element>::reverse_iterator iterL;
//		//std::vector<Element> *elementsR;
//		std::vector<Element>::iterator iterR;
//	public:
//		FracIter(): iterL(elementsL.rbegin()), iterR(elementsR.begin()) {};
//		void operator++() {
//			if (iterL != elementsL.rend())
//				iterL++;
//			else
//				iterR++;
//		};						   
//		
//		Element *current() {
//			if (iterL != elementsL.rend()) {
//				return &(*iterL);
//			}
//			else
//				return &(*iterR);
//		}
//	} fracIter;
	
	
	Breaker breaker;	//	class of the inner breaker
	Stratum *stratum;	//	pointer to the stratum which owns the fracture
	
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
	void addNewElements(const double &deltaBeta1, const double &deltaBeta2,
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
	 * near the given element 
     * @param element1 element to calculate the rotation of propagation direction
     * @return angle of rotation of propagation direction
     */
	double calcAngleOfRotation(const Element &element1) const;
};

namespace std {
	inline std::ostream& operator<<(std::ostream &os, const Fracture &frac) {
		os << "\t" << frac.getNumber()
				<< std::endl;
		return os;
	};
}

#endif	/* FRACTURE_HPP */

