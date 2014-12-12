#ifndef FLUID_HPP
#define	FLUID_HPP

#include <string>
#include "Break.hpp"

/**
 * Class for simulation of the fluid in the fracture (at this moment just by 
 * stupid setting the pressure) and for providing interaction between 
 * fluid and fracture (now just  ( pressure = - sigmaN of the element ) )
 *  
 */


class Fluid {
public:
	Fluid();
	~Fluid();
	/**
	 * Pressure of the inner fluid of the fracture is equal to
	 * - (at^2 + bt + c), t from -1 to 1, if pressureType is "polynomial",
	 * and ( -c ), if pressureType is "const" or "lag",
	 * pressureType "lag" means in tips of the fracture pressure = 0 
     * @param _a
     * @param _b
     * @param _c
     * @param _pressureType
     */
	void setType(double _a, double _b, double _c, std::string _pressureType);
	/**
	 * Calculate pressure of the fluid in the fracture and set sigmaN
	 * of elements of fracture to -pressure
     * @param breaks pointer to the middle (!) of the breaks of the fracture
     * @param NumOfCalcBrks actual number of breaks in the fracture
     */
	void calculatePressure(Break *breaks, int NumOfCalcBrks);
private:
	double a, b, c;	//	polynomial pressure parameters
	std::string pressureType;	//	type of pressure distribution along the fracture
};

#endif	/* FLUID_HPP */

