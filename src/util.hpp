#ifndef UTIL_HPP
#define	UTIL_HPP

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>

#include <cmath>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

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
void info(const T t, const Args... args) {
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
void print(const T t, const Args... args) {
	std::cerr << t << " ";
	print(args...);
}

/**
 * Print NxN System of Linear Equations A*x = b
 * @param A matrix of SLE
 * @param x sought vector in SLE
 * @param b right-hand side vector of SLE
 * @param N size of vector and matrix
 */
void printSLE(const gsl_matrix *A, const gsl_vector *x,
              const gsl_vector *b, const int N);

/**
 * Check the quality of the solution of NxN System of Linear Equations A*x = b
 * @param A matrix of SLE
 * @param x sought vector in SLE
 * @param b right-hand side vector of SLE
 * @param N size of vector and matrix
 */
void checkSLE(const gsl_matrix *A, const gsl_vector *x, 
              const gsl_vector *b, const int N);

/**
 * Calculate arctan of (a / b) taking into account the infinities
 * @param a	numerator
 * @param b	denominator
 * @return arctan (a / b)
 */
double arctan(const double &a, const double &b);

#endif	/* UTIL_HPP */