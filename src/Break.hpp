#ifndef BREAK_HPP
#define	BREAK_HPP

#include <cmath>
#include <iostream>
#include "Field.hpp"

class Break {
public:
	double Ds;	//	Break along
	double Dn;	//	Break across
	double Cx;	//	Center on x
	double Cy;	//	Center on y
	double beta;	//	Angle to x
	
	Break();
	Break(double a, double Cx, double Cy, double beta, double G, double nu);
	~Break();
	double getBs() const;
	double getBn() const;
	void setSigmaN(const double &_sigmaN);
	void calculateImpactOn(const Break &break2, double &Ass, double &Asn, double &Ans, double &Ann);
	Field calculateImpactInPoint(const double &x_glob, const double &y_glob) const;
	void setExternalImpact(Field field);
	double K1() const;
	double K2() const;
	
private:
	double G, nu;	//	Rheology parameters
	double a;	//	Half-length
	double sigmaN;	//	Impact of the inner fluid
	double sigmaS;	//	Impact of the inner fluid
	double externalSigmaN;	//	Impact of already existing fractures 
	double externalSigmaS;	//	Impact of already existing fractures
	
	double F1(const double &x, const double &y) const;
	double F2(const double &x, const double &y) const;
	double F3(const double &x, const double &y) const;
	double F4(const double &x, const double &y) const;
	double F5(const double &x, const double &y) const;
	double F6(const double &x, const double &y) const;
	double F7(const double &x, const double &y) const;
	
};

namespace std {
	inline std::ostream& operator<<(std::ostream &os, const Break &brk) {
		os << brk.Cx << "\t" << brk.Cy
				<< "\t" << brk.Ds << "\t" << brk.Dn 
				<< "\t" << brk.getBs() << "\t" << brk.getBn()  
				<< std::endl;
		return os;
	};
}

#endif	/* BREAK_HPP */

