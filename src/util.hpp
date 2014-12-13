#ifndef UTIL_HPP
#define	UTIL_HPP

#include <iostream>
#include <string>
#include <stdexcept>

void info();
template<typename T, typename... Args>
void info(T t, Args... args) {
	std::cout << t << " ";
	info(args...);
}

void print();
template<typename T, typename... Args>
void print(T t, Args... args) {
	std::cerr << t << " ";
	print(args...);
}




#endif	/* UTIL_HPP */