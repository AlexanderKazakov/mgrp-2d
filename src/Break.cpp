#include "Break.hpp"

Break::Break() {
}
	
Break::Break(double a, double Cx, double Cy, double beta, double G,
	double nu): a(a), Cx(Cx), Cy(Cy), beta(beta), G(G), nu(nu) {
	sigmaS = 0;
}

Break::~Break() {
}

void Break::setExternalImpact(Field field) {
	externalSigmaN = field.inRotatedAxis(beta).Syy;
	externalSigmaS = field.inRotatedAxis(beta).Sxy;
}

void Break::calculateImpactOf(const Break& break2, double& Ass, double& Asn,
													double& Ans, double& Ann) {
	double gamma = beta - break2.getBeta();
	double x = (Cx - break2.getCx()) * cos(break2.getBeta()) +
								(Cy - break2.getCy()) * sin(break2.getBeta());
	double y = - (Cx - break2.getCx()) * sin(break2.getBeta()) +
								(Cy - break2.getCy()) * cos(break2.getBeta());
	Ass = 2*G * ( - sin(2 * gamma) * F4(x, y) - cos(2 * gamma) * F5(x, y) - 
			y * (sin(2*gamma) * F6(x, y) - cos(2*gamma) * F7(x, y)) );
	Asn =-2*G * y * (cos(2*gamma) * F6(x, y) + sin(2*gamma) * F7(x, y));
	Ans = 2*G * ( 2 * sin(gamma) * sin(gamma) * F4(x, y) + sin(2 * gamma) * F5(x, y) - 
			y * (cos(2*gamma) * F6(x, y) + sin(2*gamma) * F7(x, y)) );
	Ann = 2*G * ( -F5(x, y) + y * (sin(2*gamma) * F6(x, y) - cos(2*gamma) * F7(x, y)) );
}

Field Break::calculateImpactInPoint(const double &x_glob, const double &y_glob) const {
	double x = (x_glob - Cx) * cos(beta) + (y_glob - Cy) * sin(beta);
	double y = (y_glob - Cy) * cos(beta) - (x_glob - Cx) * sin(beta);
	
	Field field;

	field.Ux = Ds * (-(1 - 2 * nu) * sin(beta) * F2(x, y) + 2 * (1 - nu) * cos(beta) * F3(x, y) +
			y * (sin(beta) * F4(x, y) - cos(beta) * F5(x, y))) +
			Dn * (-(1 - 2 * nu) * cos(beta) * F2(x, y) - 2 * (1 - nu) * sin(beta) * F3(x, y) -
			y * (cos(beta) * F4(x, y) + sin(beta) * F5(x, y)));
	
	field.Uy = Ds * ((1 - 2 * nu) * cos(beta) * F2(x, y) + 2 * (1 - nu) * sin(beta) * F3(x, y) -
			y * (cos(beta) * F4(x, y) + sin(beta) * F5(x, y))) +
			Dn * (-(1 - 2 * nu) * sin(beta) * F2(x, y) + 2 * (1 - nu) * cos(beta) * F3(x, y) -
			y * (sin(beta) * F4(x, y) - cos(beta) * F5(x, y)));
	
	field.Sxx = 2 * G * Ds * (2 * cos(beta) * cos(beta) * F4(x, y) + sin(2 * beta) * F5(x, y) +
			y * (cos(2 * beta) * F6(x, y) - sin(2 * beta) * F7(x, y))) +
			2 * G * Dn * (-F5(x, y) + y * (sin(2 * beta) * F6(x, y) + cos(2 * beta) * F7(x, y)));
	
	field.Syy = 2 * G * Ds * (2 * sin(beta) * sin(beta) * F4(x, y) - sin(2 * beta) * F5(x, y) -
			y * (cos(2 * beta) * F6(x, y) - sin(2 * beta) * F7(x, y))) +
			2 * G * Dn * (-F5(x, y) - y * (sin(2 * beta) * F6(x, y) + cos(2 * beta) * F7(x, y)));
	
	field.Sxy = 2 * G * Ds * (sin(2 * beta) * F4(x, y) - cos(2 * beta) * F5(x, y) +
			y * (sin(2 * beta) * F6(x, y) + cos(2 * beta) * F7(x, y))) +
			2 * G * Dn * (-y * (cos(2 * beta) * F6(x, y) - sin(2 * beta) * F7(x, y)));

	return field;
}

double Break::getBeta() const {
	return beta;
}

double Break::getCx() const {
	return Cx;
}

double Break::getCy() const {
	return Cy;
}

double Break::getDs() const {
	return Ds;
}

double Break::getDn() const {
	return Dn;
}

double Break::getBs() const {
	return sigmaS - externalSigmaS;
}

double Break::getBn() const {
	return sigmaN - externalSigmaN;
}

void Break::setDs(const double &_Ds) {
	Ds = _Ds;
}

void Break::setDn(const double &_Dn) {
	Dn = _Dn;
}

void Break::setSigmaN(const double& _sigmaN) {
	sigmaN = _sigmaN;
}

double Break::F1(const double& x, const double& y) const {
	return - ( y * (atan(y / (x - a)) - atan(y / (x + a))) - 
			(x - a) * log (sqrt( (x - a)*(x - a) + y*y)) +
			(x + a) * log (sqrt( (x + a)*(x + a) + y*y)) )
			/ ( 4 * M_PI * (1 - nu) );
}

double Break::F2(const double& x, const double& y) const {
	return ( log(sqrt((x - a)*(x - a) + y*y)) - log(sqrt((x + a)*(x + a) + y*y)) )
			/ ( 4 * M_PI * (1 - nu) );
}

double Break::F3(const double& x, const double& y) const {
	return - ( atan(y / (x - a)) - atan(y / (x + a)) ) / ( 4 * M_PI * (1 - nu) );
}

double Break::F4(const double& x, const double& y) const {
	return ( y / ((x - a)*(x - a) + y*y) - y / ((x + a)*(x + a) + y*y)  )
			/ ( 4 * M_PI * (1 - nu) );
}

double Break::F5(const double& x, const double& y) const {
	return ( (x - a) / ((x - a)*(x - a) + y*y) - (x + a) / ((x + a)*(x + a) + y*y) )
			/ ( 4 * M_PI * (1 - nu) );
}

double Break::F6(const double& x, const double& y) const {
	return ( ((x - a)*(x - a) - y*y) / ((x - a)*(x - a) + y*y) / ((x - a)*(x - a) + y*y) 
			- ((x + a)*(x + a) - y*y) / ((x + a)*(x + a) + y*y) / ((x + a)*(x + a) + y*y)  )
			/ ( 4 * M_PI * (1 - nu) );
}

double Break::F7(const double& x, const double& y) const {
	return 2 * y * ( (x - a) / ((x - a)*(x - a) + y*y) / ((x - a)*(x - a) + y*y) 
			- (x + a) / ((x + a)*(x + a) + y*y) / ((x + a)*(x + a) + y*y)  )
			/ ( 4 * M_PI * (1 - nu) );
}