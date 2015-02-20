#include "Element.hpp"

Element::Element() {
}
	
Element::Element(double a, double Cx, double Cy, 
                 double beta, double G, double nu):
                 a(a), Cx(Cx), Cy(Cy),
                 beta(beta), G(G), nu(nu), 
                 sigmaS(0) {
}

Element::~Element() {
}

double Element::getBs() const {
	return sigmaS - externalSigmaS;
}

double Element::getBn() const {
	return sigmaN - externalSigmaN;
}

double Element::getA() const {
	return a;
}

void Element::setSigmaN(const double& _sigmaN) {
	sigmaN = _sigmaN;
}

void Element::setExternalImpact(const Field &field) {
	externalSigmaN = field.inRotatedAxis(beta).Syy;
	externalSigmaS = field.inRotatedAxis(beta).Sxy;
}

void Element::calculateImpactOn(const Element& lmnt2,
                                double& Ass, double& Asn,
                                double& Ans, double& Ann) const {
	double gamma = lmnt2.beta - beta;
	// x and y are coordinates of lmnt2 in the system 
	// of coordinates of this element.
	double x = (lmnt2.Cx - Cx) * cos(beta) + (lmnt2.Cy - Cy) * sin(beta);
	double y = - (lmnt2.Cx - Cx) * sin(beta) + (lmnt2.Cy - Cy) * cos(beta);
	
	Ass = 2 * G * ( - sin(2 * gamma) * F4(x, y) 
	                - cos(2 * gamma) * F5(x, y) 
	                - y * (sin(2*gamma) * F6(x, y) 
	                - cos(2*gamma) * F7(x, y)) );
	
	Asn = - 2 * G * y * (  cos(2*gamma) * F6(x, y) 
	                     + sin(2*gamma) * F7(x, y));
	
	Ans = 2 * G * ( 2 * sin(gamma) * sin(gamma) * F4(x, y) 
	                + sin(2 * gamma) * F5(x, y) 
	                - y * (cos(2*gamma) * F6(x, y) 
	                + sin(2*gamma) * F7(x, y)) );
	
	Ann = 2 * G * ( - F5(x, y) + y * (sin(2*gamma) * F6(x, y)
			                   - cos(2*gamma) * F7(x, y)) );
}

Field Element::calculateImpactInPoint(const double &x_glob, 
	                                  const double &y_glob) const {
	// x and y are coordinates of given point in the system 
	// of coordinates of this element.
	double x = (x_glob - Cx) * cos(beta) + (y_glob - Cy) * sin(beta);
	double y = (y_glob - Cy) * cos(beta) - (x_glob - Cx) * sin(beta);
	
	Field field;

	field.Ux = Ds * (-(1 - 2 * nu) * sin(beta) * F2(x, y) + 
	                 2 * (1 - nu) * cos(beta) * F3(x, y) +
	                 y * (sin(beta) * F4(x, y) - cos(beta) * F5(x, y)))
	         + Dn * (-(1 - 2 * nu) * cos(beta) * F2(x, y) - 
	                 2 * (1 - nu) * sin(beta) * F3(x, y) -
	                 y * (cos(beta) * F4(x, y) + sin(beta) * F5(x, y)));
	
	field.Uy = Ds * ((1 - 2 * nu) * cos(beta) * F2(x, y) + 
	                 2 * (1 - nu) * sin(beta) * F3(x, y) -
	                 y * (cos(beta) * F4(x, y) + sin(beta) * F5(x, y)))
	         + Dn * (-(1 - 2 * nu) * sin(beta) * F2(x, y) + 
	                 2 * (1 - nu) * cos(beta) * F3(x, y) -
	                 y * (sin(beta) * F4(x, y) - cos(beta) * F5(x, y)));
	
	field.Sxx = 2 * G * Ds * ( 2 * cos(beta) * cos(beta) * F4(x, y) 
	                           + sin(2 * beta) * F5(x, y) 
	                           + y * (cos(2 * beta) * F6(x, y) 
	                           - sin(2 * beta) * F7(x, y))) 
	          + 2 * G * Dn * ( - F5(x, y) + y * (sin(2 * beta) * F6(x, y) 
	                                      + cos(2 * beta) * F7(x, y)));
	
	field.Syy = 2 * G * Ds * ( 2 * sin(beta) * sin(beta) * F4(x, y) 
	                           - sin(2 * beta) * F5(x, y) 
	                           - y * (cos(2 * beta) * F6(x, y) 
	                           - sin(2 * beta) * F7(x, y))) 
	          + 2 * G * Dn * ( - F5(x, y) - y * (sin(2 * beta) * F6(x, y) 
	                                      + cos(2 * beta) * F7(x, y)));
	
	field.Sxy = 2 * G * Ds * (  sin(2 * beta) * F4(x, y) 
	                          - cos(2 * beta) * F5(x, y) 
	                          + y * (sin(2 * beta) * F6(x, y) 
	                          + cos(2 * beta) * F7(x, y))) 
	          + 2 * G * Dn * ( - y * (cos(2 * beta) * F6(x, y) 
			                   - sin(2 * beta) * F7(x, y)));

	return field;
}

double Element::K1() const {
	return G / 4 / (1 - nu) * sqrt(2 * M_PI / a) * Dn;
}

double Element::K2() const {
	return G / 4 / (1 - nu) * sqrt(2 * M_PI / a) * Ds;
}

double Element::F1(const double& x, const double& y) const {
	return - ( y * (atan(y / (x - a)) - atan(y / (x + a))) - 
	           (x - a) * log (sqrt( (x - a)*(x - a) + y*y)) +
	           (x + a) * log (sqrt( (x + a)*(x + a) + y*y)) )
	         / ( 4 * M_PI * (1 - nu) );
}

double Element::F2(const double& x, const double& y) const {
	return (   log(sqrt((x - a)*(x - a) + y*y)) 
	         - log(sqrt((x + a)*(x + a) + y*y)) )
	       / ( 4 * M_PI * (1 - nu) );
}

double Element::F3(const double& x, const double& y) const {
	return - ( atan(y / (x - a)) - atan(y / (x + a)) ) 
	         / ( 4 * M_PI * (1 - nu) );
}

double Element::F4(const double& x, const double& y) const {
	return ( y / ((x - a)*(x - a) + y*y) - y / ((x + a)*(x + a) + y*y)  )
	       / ( 4 * M_PI * (1 - nu) );
}

double Element::F5(const double& x, const double& y) const {
	return (   (x - a) / ((x - a)*(x - a) + y*y) 
	         - (x + a) / ((x + a)*(x + a) + y*y) )
	       / ( 4 * M_PI * (1 - nu) );
}

double Element::F6(const double& x, const double& y) const {
	return (  ((x - a)*(x - a) - y*y)
	          / ((x - a)*(x - a) + y*y) / ((x - a)*(x - a) + y*y) 
			- ((x + a)*(x + a) - y*y)
	          / ((x + a)*(x + a) + y*y) / ((x + a)*(x + a) + y*y)  )
	       / ( 4 * M_PI * (1 - nu) );
}

double Element::F7(const double& x, const double& y) const {
	return 2 * y * (  (x - a) 
	                  / ((x - a)*(x - a) + y*y) / ((x - a)*(x - a) + y*y) 
	                - (x + a) 
	                  / ((x + a)*(x + a) + y*y) / ((x + a)*(x + a) + y*y)  )
	               / ( 4 * M_PI * (1 - nu) );
}