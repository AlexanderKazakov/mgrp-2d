#ifndef STRATUM_HPP
#define	STRATUM_HPP

#include <vector>
#include <mgl2/mgl.h>
#include "Fracture.hpp"
#include "Field.hpp"
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
	 * Pressure of the inner fluid of the fracture is equal to
	 * - (at^2 + bt + c), t from -1 to 1, if pressureType is polynomial,
	 * and ( -c ), if pressureType is const
     * @param number index number of fracture
     * @param x	initial position of fracture
     * @param y initial position of fracture
     * @param beta initial angle of fracture to x-axis
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
	void addFracture(int number, double x, double y, double beta, double h_length,
		int numOfBreaks, double a, double b, double c, std::string pressureType,
		std::string tip, std::string rotation);
	/**
	 * Set rheology parameters of stratum
     * @param _G shear modulus
     * @param _nu Poisson's ratio
     */
	void setRheology(double _G, double _nu);
	/**
	 * Get rheology parameters from stratum
     * @param _G shear modulus
     * @param _nu Poisson's ratio
     */
	void getRheology(double &_G, double &_nu);
	/**
	 * Set external ( "at infinity" ) stresses
     * @param _Sxx
     * @param _Sxy
     * @param _Syy
     */
	void setStresses(double _Sxx, double _Sxy, double _Syy);
	/**
	 * Set ranges for visualization graph
     * @param _Xmin
     * @param _Xmax
     * @param _Ymin
     * @param _Ymax
     */
	void setRanges(double _Xmin, double _Xmax, double _Ymin, double _Ymax);
	/**
	 * Reserve enough memory for numberOfFractures fractures.
	 * It's very important for correct work of the program!
     * @param numberOfFractures 
     */
	void reserve(int numberOfFractures);
	/**
	 * Calculate the next one fracture supposing that all previous 
	 * fractures stay constant
     */
	void calculate();
	/**
	 * Calculate field in point (x, y) caused by existing fractures and
	 * external stresses
     * @param x	x-coord of point
     * @param y y-coord of point
     * @return Actual Field in (x, y)
     */
	Field calculateImpactInPoint(const double &x, const double &y);
	/**
	 * Create enterprise graphics
     */
	void visualize(); 
	/**
	 * Draw displacement discontinuities for one fracture
     */
	void drawDisplacements();
private:
	double G;	//	Shear modulus
	double nu;	//	Poisson's ratio
	double E;	//	Young's modulus	E = 2 * G * (1 + nu)
	double Xmin, Xmax, Ymin, Ymax;	//	Ranges for plot
	double Sxx;	//	
	double Sxy;	//	External stresses
	double Syy;	//
	std::vector<Fracture> fractures;
	std::vector<Fracture>::iterator currentFracture;
						// the iterator on the fracture that is calculated now
	/**
	 * Draw the actual field
     * @param gr graph to draw on
     */
	void drawField(mglGraph &gr);
	/**
	 * Draw actual fractures
     * @param gr graph to draw on
     */
	void drawFractures(mglGraph &gr);
	/**
	 * Draw main stress' directions
     * @param gr graph to draw on
     */
	void drawStressDirections(mglGraph &gr);
};

#endif	/* STRATUM_HPP */

