#include "Fracture.hpp"

#include <c++/4.9/vector>

Fracture::Fracture() {
}

Fracture::Fracture(Stratum *stratum, int number, int numOfLmntsL, int numOfLmntsR,
                   double halfLengthOfLmnts, double _a, double _b, double _c,
                   std::string pressureType, std::string tip,
                   std::string rotation): 
                   stratum(stratum), number(number), numOfLmntsL(numOfLmntsL),
                   numOfLmntsR(numOfLmntsR), a(halfLengthOfLmnts),
                   rotation(rotation), tip(tip) {
	fractionIsStoppedL = false;
	fractionIsStoppedR = false;
	double _G, _nu;
	stratum->getRheology(_G, _nu);
	G = _G;	nu = _nu;
	breaker.setType(_a, _b, _c, pressureType);
}

Fracture::~Fracture() {
}

void Fracture::allocateLmnts(double x, double y, double beta) {
	lmntsL.reserve(numOfLmntsL);
	lmntsR.reserve(numOfLmntsR);
	
	// Add three small tip elements to the corners
	lmntsR.push_back(Element(5 * a / 9, x + 5 * a / 9 * cos(beta), 
	                               y + 5 * a / 9 * sin(beta), beta, G, nu));
	lmntsR.push_back(Element(a / 3, x + 13 * a / 9 * cos(beta), 
	                           y + 13 * a / 9 * sin(beta), beta, G, nu));
	lmntsR.push_back(Element(a / 9, x + 17 * a / 9 * cos(beta), 
	                           y + 17 * a / 9 * sin(beta), beta, G, nu));
	lmntsL.push_back(Element(5 * a / 9, x - 5 * a / 9 * cos(beta), 
	                               y - 5 * a / 9 * sin(beta), beta, G, nu));
	lmntsL.push_back(Element(a / 3, x - 13 * a / 9 * cos(beta), 
	                           y - 13 * a / 9 * sin(beta), beta, G, nu));
	lmntsL.push_back(Element(a / 9, x - 17 * a / 9 * cos(beta), 
	                           y - 17 * a / 9 * sin(beta), beta, G, nu));
	
	halfLengthL = halfLengthR = 2*a;
}

Field Fracture::calculateImpactInPoint(const double &x, const double &y) const {
	Field field;
	for (auto lmnt = lmntsL.begin(); lmnt != lmntsL.end(); lmnt++) {
		field += lmnt->calculateImpactInPoint(x, y);
	}
	for (auto lmnt = lmntsR.begin(); lmnt != lmntsR.end(); lmnt++) {
		field += lmnt->calculateImpactInPoint(x, y);
	}
	return field;
}

void Fracture::calculate() {
	info("Starting calculation of fracture number", number, "...");
	for (auto lmnt = lmntsL.begin(); lmnt != lmntsL.end(); lmnt++) 
		lmnt->setExternalImpact
	                     (stratum->calculateImpactInPoint(lmnt->Cx, lmnt->Cy));
	for (auto lmnt = lmntsR.begin(); lmnt != lmntsR.end(); lmnt++) 
		lmnt->setExternalImpact
	                     (stratum->calculateImpactInPoint(lmnt->Cx, lmnt->Cy));
	breaker.calculatePressure(lmntsL, lmntsR);
	
	calculateLmnts();
	if (fractionIsStoppedL * fractionIsStoppedR)
		info("Fracture number", number, "is stopped.");
	else if (fractionIsStoppedL)
		info("Fracture number", number, "is stopped at the left corner.");
	else if (fractionIsStoppedR)
		info("Fracture number", number, "is stopped at the right corner.");
				
	while ( ! (fractionIsStoppedL * fractionIsStoppedR) ) {	
		double deltaBeta1 = calcAngleOfRotation(lmntsL.back());
		double deltaBeta2 = calcAngleOfRotation(lmntsR.back());
		
		// TODO - add the opportunity to work with a usual const tip element
		grow(deltaBeta1, deltaBeta2);
		if (rotation == "predictor-corrector") {
			//	Clarifying the direction of fracture's growth like it is
			//	usually done in "predictor-corrector" methods
			calculateLmnts();
			deltaBeta1 = (calcAngleOfRotation(lmntsL.back()) + deltaBeta1) / 2;
			deltaBeta2 = (calcAngleOfRotation(lmntsR.back()) + deltaBeta2) / 2;		
			replaceTipElements(deltaBeta1, deltaBeta2);
		}
		calculateLmnts();
		if (fractionIsStoppedL * fractionIsStoppedR)
			info("Fracture number", number, "is stopped.");
		else if (fractionIsStoppedL)
			info("Fracture number", number, "is stopped at the left corner.");
		else if (fractionIsStoppedR)
			info("Fracture number", number, "is stopped at the right corner.");
		
		if ( (lmntsL.size() >= numOfLmntsL) && (lmntsR.size() >= numOfLmntsR) )
			break;
	}
	
	info("Fracture number", number, "is calculated.");
}

void Fracture::calculateLmnts() {
	double Ass, Asn, Ans, Ann;
	int N = 2 * (lmntsL.size() + lmntsR.size());
	gsl_matrix *A = gsl_matrix_alloc(N, N);
	gsl_vector *b = gsl_vector_alloc(N);
	gsl_vector *x = gsl_vector_alloc(N);

	// Filling in the matrix and right-hand-side vector of system
	// of linear equations on displacement discontinuities
	int i = 0;
	for (auto lmnt1 = lmntsL.rbegin(); lmnt1 != lmntsL.rend(); lmnt1++) {
		int j = 0;
		for (auto lmnt2 = lmntsL.rbegin(); lmnt2 != lmntsL.rend(); lmnt2++) {
			lmnt2->calculateImpactOn(*lmnt1, Ass, Asn, Ans, Ann);
			gsl_matrix_set(A, i, j, Ass);
			gsl_matrix_set(A, i, j + 1, Asn);
			gsl_matrix_set(A, i + 1, j, Ans);
			gsl_matrix_set(A, i + 1, j + 1, Ann);
			j += 2;
		}
		for (auto lmnt2 = lmntsR.begin(); lmnt2 != lmntsR.end(); lmnt2++) {
			lmnt2->calculateImpactOn(*lmnt1, Ass, Asn, Ans, Ann);
			gsl_matrix_set(A, i, j, Ass);
			gsl_matrix_set(A, i, j + 1, Asn);
			gsl_matrix_set(A, i + 1, j, Ans);
			gsl_matrix_set(A, i + 1, j + 1, Ann);
			j += 2;
		}
		gsl_vector_set(b, i, lmnt1->getBs());
		gsl_vector_set(b, i + 1, lmnt1->getBn());
		i += 2;
	}
	for (auto lmnt1 = lmntsR.begin(); lmnt1 != lmntsR.end(); lmnt1++) {
		int j = 0;
		for (auto lmnt2 = lmntsL.rbegin(); lmnt2 != lmntsL.rend(); lmnt2++) {
			lmnt2->calculateImpactOn(*lmnt1, Ass, Asn, Ans, Ann);
			gsl_matrix_set(A, i, j, Ass);
			gsl_matrix_set(A, i, j + 1, Asn);
			gsl_matrix_set(A, i + 1, j, Ans);
			gsl_matrix_set(A, i + 1, j + 1, Ann);
			j += 2;
		}
		for (auto lmnt2 = lmntsR.begin(); lmnt2 != lmntsR.end(); lmnt2++) {
			lmnt2->calculateImpactOn(*lmnt1, Ass, Asn, Ans, Ann);
			gsl_matrix_set(A, i, j, Ass);
			gsl_matrix_set(A, i, j + 1, Asn);
			gsl_matrix_set(A, i + 1, j, Ans);
			gsl_matrix_set(A, i + 1, j + 1, Ann);
			j += 2;
		}
		gsl_vector_set(b, i, lmnt1->getBs());
		gsl_vector_set(b, i + 1, lmnt1->getBn());
		i += 2;
	}
	
	// Solving the SLE by gsl-library
	gsl_permutation *p = gsl_permutation_alloc(N);
	int signum;		
	gsl_linalg_LU_decomp(A, p, &signum);
	gsl_linalg_LU_solve(A, p, b, x);

//	for (int p = 0; p < N; p++) {
//		for (int q = 0; q < N; q++) {
//			std::cout << gsl_matrix_get(A, p, q) << "\t";
//		}
//		std::cout << "\t*" << gsl_vector_get(x, p)  << "\t=" << gsl_vector_get(b, p) << std::endl;
//	}
//	std::cout << std::endl;
	
	i = 0;
	for (auto lmnt1 = lmntsL.rbegin(); lmnt1 != lmntsL.rend(); lmnt1++) {
		lmnt1->Ds = gsl_vector_get(x, i);
		lmnt1->Dn = gsl_vector_get(x, i + 1);
		i += 2;
	}
	for (auto lmnt1 = lmntsR.begin(); lmnt1 != lmntsR.end(); lmnt1++) {
		lmnt1->Ds = gsl_vector_get(x, i);
		lmnt1->Dn = gsl_vector_get(x, i + 1);
		i += 2;
	}
	gsl_matrix_free(A);
	gsl_vector_free(x);
	gsl_vector_free(b);

	if (lmntsL.back().Dn > 0)
		fractionIsStoppedL = true;
	if (lmntsR.back().Dn > 0)
		fractionIsStoppedR = true;
	
}

double Fracture::calcAngleOfRotation(const Element &lmnt1) const {
	double K1 = - lmnt1.Dn;
	double K2 = - lmnt1.Ds;
	Field tmp;
	double beta = 2 * tmp.arctan(- 2 * K2, K1 + sqrt(K1 * K1 + 8 * K2 * K2) );
	return beta;
}

void Fracture::grow(const double &deltaBeta1, const double &deltaBeta2) {
	if ( (lmntsL.size() == 3) && (lmntsR.size() == 3) ) {
		double beta0 = lmntsL.back().beta;
		double x0 = lmntsL.front().Cx + 5 * a / 9 * cos(beta0);
		double y0 = lmntsL.front().Cy + 5 * a / 9 * sin(beta0);
		lmntsL.pop_back(); lmntsL.pop_back(); lmntsL.pop_back();
		lmntsR.pop_back(); lmntsR.pop_back(); lmntsR.pop_back();
		lmntsL.push_back(Element(a, x0 - a * cos(beta0), 
	                                y0 - a * sin(beta0), beta0, G, nu));
		lmntsL.back().setExternalImpact(stratum->calculateImpactInPoint
	                                     (lmntsL.back().Cx, lmntsL.back().Cy));
		lmntsR.push_back(Element(a, x0 + a * cos(beta0), 
	                                y0 + a * sin(beta0), beta0, G, nu));
		lmntsR.back().setExternalImpact(stratum->calculateImpactInPoint
	                                     (lmntsL.back().Cx, lmntsL.back().Cy));
		addNewLmnts(deltaBeta1, deltaBeta2, 5 * a / 9);
		addNewLmnts(0, 0, a / 3);
		addNewLmnts(0, 0, a / 9);
	} else {
		double beta1 = lmntsL.back().beta;
		double beta2 = lmntsR.back().beta;
		lmntsL.pop_back(); lmntsL.pop_back(); lmntsL.pop_back();
		lmntsR.pop_back(); lmntsR.pop_back(); lmntsR.pop_back();
		halfLengthL -= 2*a; halfLengthR -= 2*a;
		double oldDeltaBeta1 = beta1 - lmntsL.back().beta;
		double oldDeltaBeta2 = beta2 - lmntsR.back().beta;
		
		addNewLmnts(oldDeltaBeta1, oldDeltaBeta2, a);
		addNewLmnts(deltaBeta1, deltaBeta2, 5 * a / 9);
		addNewLmnts(0, 0, a / 3);
		addNewLmnts(0, 0, a / 9);
	}
	
	breaker.calculatePressure(lmntsL, lmntsR);
}

void Fracture::addNewLmnts(const double& deltaBeta1, const double& deltaBeta2, 
                                                     const double& half_length) {
	double beta1 = lmntsL.back().beta + deltaBeta1;
	double x1 = lmntsL.back().Cx - lmntsL.back().getA() * cos(lmntsL.back().beta) 
	                             - half_length * cos(beta1);
	double y1 = lmntsL.back().Cy - lmntsL.back().getA() * sin(lmntsL.back().beta)
	                             - half_length * sin(beta1);

	lmntsL.push_back(Element(half_length, x1, y1, beta1, G, nu));
	lmntsL.back().setExternalImpact(stratum->calculateImpactInPoint
	                                         (lmntsL.back().Cx, lmntsL.back().Cy));

	beta1 = lmntsR.back().beta + deltaBeta2;
	x1 = lmntsR.back().Cx + lmntsR.back().getA() * cos(lmntsR.back().beta) 
	                      + half_length * cos(beta1);
	y1 = lmntsR.back().Cy + lmntsR.back().getA() * sin(lmntsR.back().beta)
	                      + half_length * sin(beta1);

	lmntsR.push_back(Element(half_length, x1, y1, beta1, G, nu));
	lmntsR.back().setExternalImpact(stratum->calculateImpactInPoint
	                                         (lmntsR.back().Cx, lmntsR.back().Cy));
	
	halfLengthL += 2 * half_length;
	halfLengthR += 2 * half_length;
}

void Fracture::replaceTipElements(const double& deltaBeta1, const double& deltaBeta2) {
	lmntsL.pop_back(); lmntsL.pop_back(); lmntsL.pop_back();
	lmntsR.pop_back(); lmntsR.pop_back(); lmntsR.pop_back();
	halfLengthL -= 2*a; halfLengthR -= 2*a;
	addNewLmnts(deltaBeta1, deltaBeta2, 5 * a / 9);
	addNewLmnts(0, 0, a / 3);
	addNewLmnts(0, 0, a / 9);
	breaker.calculatePressure(lmntsL, lmntsR);
}

int Fracture::getNumOfLmntsL() const {
	return lmntsL.size();
}

int Fracture::getNumOfLmntsR() const {
	return lmntsR.size();
}

void Fracture::getPointsForPlot(double* x, double* y) const {
	x[0] = lmntsL.back().Cx - lmntsL.back().getA() * cos(lmntsL.back().beta);
	y[0] = lmntsL.back().Cy - lmntsL.back().getA() * sin(lmntsL.back().beta);

	int i = 1;
	for (auto lmnt = lmntsL.rbegin(); lmnt != lmntsL.rend(); lmnt++) {
		x[i] = lmnt->Cx + lmnt->getA() * cos(lmnt->beta);
		y[i] = lmnt->Cy + lmnt->getA() * sin(lmnt->beta);
		i++;
	}
	for (auto lmnt = lmntsR.begin(); lmnt != lmntsR.end(); lmnt++) {
		x[i] = lmnt->Cx + lmnt->getA() * cos(lmnt->beta);
		y[i] = lmnt->Cy + lmnt->getA() * sin(lmnt->beta);
		i++;
	}
}

void Fracture::getPointsForDisplacementPlot(double* x, double* v) const {
	int i = 0;
	for (auto lmnt = lmntsL.rbegin(); lmnt != lmntsL.rend(); lmnt++) {
		x[i] = lmnt->Cx;
		v[i] = lmnt->Dn / 2;
		i++;
	}
	for (auto lmnt = lmntsR.begin(); lmnt != lmntsR.end(); lmnt++) {
		x[i] = lmnt->Cx;
		v[i] = lmnt->Dn / 2;
		i++;
	}
}

double Fracture::getHalfLengthL() const {
	return halfLengthL;
}

double Fracture::getHalfLengthR() const {
	return halfLengthR;
}
