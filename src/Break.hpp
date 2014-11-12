#ifndef BREAK_HPP
#define	BREAK_HPP

#include <cmath>
#include <iostream>
#include "Field.hpp"

class Break {
public:
	Break();
	Break(int number, double a, double Cx, double Cy, double beta, double G, double nu,
										double sigmaN, double sigmaS);
	~Break();
	int getNumber() const;
	double getBeta() const;
	double getCx() const;
	double getCy() const;
	double getBs() const;
	double getBn() const;
	double getDs() const;
	double getDn() const;
	void setDs(const double &_Ds);
	void setDn(const double &_Dn);
	bool operator==(const Break &other) const;
	bool operator!=(const Break &other) const;
	bool operator<(const Break &other) const;
	void calculateImpactOf(const Break &break2, double &Ass, double &Asn, double &Ans, double &Ann);
	Field calculateImpactInPoint(const double &x_glob, const double &y_glob) const;
	void setExternalImpact(Field field);
	
private:
	double G, nu;	//	Rheology parameters
	double Ds;	//	Break along
	double Dn;	//	Break across
	double a;	//	Half-length
	double Cx;	//	Center on x
	double Cy;	//	Center on y
	double beta;	//	Angle to x
	double sigmaN;
	double sigmaS;
	double externalSigmaN;	//	Impact of already existing fractures 
	double externalSigmaS;	//	Impact of already existing fractures
	int number;
	
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
		os //<< brk.getCx() << "\t" << brk.getCy() 
				<< "\t" << brk.getDs() 
				<< "\t" << brk.getDn() //<< std::endl << brk.getBs()
				//<< "\t" << brk.getBn()  
				<< std::endl;
		return os;
	};
}

#endif	/* BREAK_HPP */

