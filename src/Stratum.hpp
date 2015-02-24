#ifndef STRATUM_HPP
#define	STRATUM_HPP

#include <vector>

#include <mgl2/mgl.h>

#include "Fracture.hpp"
#include "Field.hpp"
#include "Breaker.hpp"
#include "util.hpp"

class Fracture;

/**
 * Class that contains all conditions and fractures and manipulates
 * all the actions including launching of fractures' computation and
 * visualization of the results.
 * 
 */


class Stratum {
public:
	Stratum();
	~Stratum();
	/**
	 * Add to the subsequent computation another fracture.
	 * Pressure of the inner breaker of the fracture is equal to
	 * - (at^2 + bt + c), t from -1 to 1, if pressureType is "polynomial",
	 * ( -c ), if pressureType is "const"
	 * ( -c ), t from -a to a, 0 elsewhere, 
	 * if pressureType is "lag" (means in tips of the fracture pressure = 0) 
     * @param number index number of fracture
     * @param x	initial position of fracture
     * @param y initial position of fracture
     * @param beta initial angle of fracture to x-axis
     * @param h_length half-length of boundary elements of fracture
     * @param numOfElements required number of boundary elements at the end of
	 * computation
     * @param a pressure parameter
     * @param b pressure parameter
     * @param c pressure parameter
     * @param pressureType type of pressure of the inner breaker
     * @param tip rule to use special boundary element at the tip of the
	 * fracture, "const" or TODO
     */
	void addFracture(int number, int numOfElements, 
                     double x, double y, double beta, double halfLengthOfElements,
                     double a, double b, double c, 
                     std::string pressureType, std::string tip);
	/**
	 * Set rheology parameters of stratum
     * @param _G shear modulus
     * @param _nu Poisson's ratio
     */
	void setRheology(const double &_G, const double &_nu);
	/**
	 * Get rheology parameters from stratum
     * @param _G shear modulus
     * @param _nu Poisson's ratio
     */
	void getRheology(double &_G, double &_nu) const;
	/**
	 * Set external ( "at infinity" ) stresses
     * @param _Sxx
     * @param _Sxy
     * @param _Syy
     */
	void setStresses(const double &_Sxx, const double &_Sxy, const double &_Syy);
	/**
	 * Set ranges for visualization graph
     * @param _Xmin
     * @param _Xmax
     * @param _Ymin
     * @param _Ymax
     */
	void setRanges(const double &_Xmin, const double &_Xmax, 
	               const double &_Ymin, const double &_Ymax);
	/**
	 * Set sequence of the calculation of fractures
     * @param _sequence "series" - one by one without feedback from next 
	 * to previous, "parallel" - all together, "series with feedback" - 
	 * one by one with feedback from next to previous
	 */
	void setSequence(const std::string _sequence);
	/**
	 * @param rotation rule to use splitting the process of fracture growing 
	 * on actual configuration's calculation and calculation of direction 
	 * of the fractures' growth, "predictor" or "predictor-corrector"
     */
	void setRotation(const std::string _rotation);
	/**
	 * Reserve enough memory for numberOfFractures fractures in 
	 * vector fractures
     * @param numberOfFractures 
     */
	void reserve(const int numberOfFractures);
	/**
	 * Calculate all required in task
     */
	void calculateTask();
	/**
	 * Calculate field in point (x, y) caused by already calculated ( from 
	 * fractures.begin() to (--beginFracture) ) fractures 
	 * and external (at infinity) stresses in the stratum
     * @param x	x-coord of point
     * @param y y-coord of point
     * @return Actual Field in (x, y)
     */
	Field calculateImpactInPoint(const double &x, const double &y) const;
	/**
	 * Create enterprise graphics
     */
	void visualize() const; 
	/**
	 * Draw displacement discontinuities for one fracture for
	 * testing purposes
     */
	void drawDisplacements() const;
private:
	double G; // Shear modulus
	double nu; // Poisson's ratio
	double Xmin, Xmax, Ymin, Ymax; // Ranges for plot
	// External stresses
	double Sxx;
	double Sxy; 
	double Syy;
	std::vector<Fracture> fractures; // all fractures
	// the iterator on the first fracture of those which are calculated currently
	std::vector<Fracture>::iterator beginFracture;
	// the iterator on the last fracture of those which are calculated currently
	std::vector<Fracture>::iterator endFracture;
	// for visualization purpose, in test2.py
	Breaker breakerOfFirstFracture;
	// rule on the sequence of calculation of several fractures,
	// "series", "parallel" or "series with feedback"
	std::string sequence;
	// rule to use splitting on actual configuration's calculation 
	// and calculation of direction of the fracture's growth
	std::string rotation;
	/**
	 * Calculate the fractures from beginFracture to endFracture 
	 * those grow together at the same time
     */
	void calculate();
	/**
	 * Fill in the system of linear equations on displacement 
	 * discontinuities of elements (for all fractures 
	 * from beginFracture to endFracture)
	 * A * x = b, 
	 * solve it by GSL library and set new values
	 * of displacement discontinuity to elements.
	 * 
	 * It's the main idea of the program. It uses displacement discontinuity 
	 * boundary method proposed by Crouch.
     */
	void calculateElements();
	/**
	 * Draw actual fractures
     * @param gr graph to draw on
     */
	void drawFractures(mglGraph &gr) const;
	/**
	 * Draw the actual field
     * @param gr graph to draw on
     */
	void drawField(mglGraph &gr) const;
	/**
	 * Draw main stress' directions
     * @param gr graph to draw on
     */
	void drawStressDirections(mglGraph &gr) const;
};

#endif	/* STRATUM_HPP */

