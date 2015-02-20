#ifndef UTIL_HPP
#define	UTIL_HPP

#include <iostream>
#include <string>
#include <stdexcept>

/**
 * Function for logging, bottom of recursion
 */
void info();

/**
 * Function for logging
 * @param t thing to print
 * @param args things to print
 */
template<typename T, typename... Args>
void info(T t, Args... args) {
	std::cout << t << " ";
	info(args...);
}

/**
 * Function for debugging, bottom of recursion
 */
void print();
/**
 * Function for debugging
 * @param t thing to print
 * @param args things to print
 */
template<typename T, typename... Args>
void print(T t, Args... args) {
	std::cerr << t << " ";
	print(args...);
}




#endif	/* UTIL_HPP */