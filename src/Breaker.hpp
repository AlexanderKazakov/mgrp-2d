#ifndef BREAKER_HPP
#define	BREAKER_HPP

#include <string>
#include <vector>

#include "Element.hpp"
#include "Fracture.hpp"
#include "util.hpp"


class Fracture;


/**
 * Class for simulation of the breaker liquid in the fracture 
 * (at this moment just by stupid setting the pressure)
 * and for providing interaction between 
 * breaker and fracture (now just  ( pressure = - sigmaN of the element ) )
 *  
 */


class Breaker {
public:
	Breaker();
	~Breaker();
	/**
	 * Pressure of the inner breaker of the fracture is equal to
	 * - (at^2 + bt + c), t from -1 to 1, if pressureType is "polynomial",
	 * ( -c ), if pressureType is "const"
	 * ( -c ), t from -a to a, 0 elsewhere, 
	 * if pressureType is "lag" (means in tips of the fracture pressure = 0) 
     * @param _a
     * @param _b
     * @param _c
     * @param _pressureType
     */
	void setType(double _a, double _b, double _c, std::string _pressureType);
	/**
	 * Get parameters from Breaker
	 * Pressure of the inner breaker of the fracture is equal to
	 * - (at^2 + bt + c), t from -1 to 1, if pressureType is "polynomial",
	 * ( -c ), if pressureType is "const"
	 * ( -c ), t from -a to a, 0 elsewhere, 
	 * if pressureType is "lag" (means in tips of the fracture pressure = 0) 
     * @param _a 
     * @param _b
     * @param _c
     * @param _pressureType
     */
	void getType(double &_a, double &_b, double &_c, 
                 std::string &_pressureType) const;
	/**
	 * Calculate pressure of the breaker in the fracture and set sigmaN
	 * of elements of fracture to -pressure
     * @param elements pointer to the middle (!) of the elements of the fracture
     * @param numOfCalcElements actual number of elements in the fracture
     */
	void calculatePressure(Fracture *frac);
private:
	// indexes to mark points of fracture at which the injection of breaker
	// was stopped
	int leftN, rightN;
	// End of Injection Multiplier. Used to take into account decrease in 
	// pressure at the end of injection
	const double eim = 1.0;
	double a, b, c; // pressure parameters
	std::string pressureType; // type of pressure distribution along the fracture
};

#endif	/* BREAKER_HPP */

