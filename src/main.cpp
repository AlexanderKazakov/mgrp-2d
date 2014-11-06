#include <cstdlib>
#include <iostream>
#include "Engine.hpp"

int main(int argc, char** argv) {
	
	Engine engine;
	
	//  Reading parameters from file
	engine.loadTask();

	//  Starting calculations
	engine.calculate();

	//	Making plots
	engine.visualize();
	
	return 0;
}
