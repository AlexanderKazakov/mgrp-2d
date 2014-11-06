#ifndef ENGINE_HPP
#define	ENGINE_HPP

#include <iostream>
#include <string>
#include <tinyxml.h>
#include"Stratum.hpp"
#include"Fracture.hpp"

class Engine {
public:
    Engine();
    ~Engine();
	void loadTask();
	void calculate();
	void visualize();
private:
	Stratum stratum;

};

#endif	/* ENGINE_HPP */

