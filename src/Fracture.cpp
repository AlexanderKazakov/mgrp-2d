#include "Fracture.hpp"

Fracture::Fracture() {
	lmnts = NULL;
}

Fracture::Fracture(Stratum *stratum, int number, int numOfLmnts,
                   double halfLengthOfLmnts, double _a, double _b, double _c,
                   std::string pressureType, std::string tip,
                   std::string rotation): 
                   stratum(stratum), number(number), numOfLmnts(numOfLmnts),
                   a(halfLengthOfLmnts),
                   rotation(rotation), tip(tip) {
	lmnts = NULL;
	fractionIsStopped = false;
	double _G, _nu;
	stratum->getRheology(_G, _nu);
	G = _G;	nu = _nu;
	front = back = middle = numOfLmnts/2;
	breaker.setType(_a, _b, _c, pressureType);
	numOfCalcLmnts = 1;
}

Fracture::~Fracture() {
	if (lmnts) {
		// Deletion of dynamical allocated memory (!)
		delete [] lmnts;
	}
}

void Fracture::allocateLmnts(double x, double y, double beta) {
	// Dynamical memory allocation (!)
	lmnts = new Element[numOfLmnts];
	lmnts[middle] = Element(a, x, y, beta, G, nu);
	halfLength = a;
}

Field Fracture::calculateImpactInPoint(const double &x, const double &y) const {
	Field field;
	for (int i = front; i <= back; i++)
		field += lmnts[i].calculateImpactInPoint(x, y);
	return field;
}

void Fracture::calculate() {
	info("Starting calculation of fracture number", number, "...");
	lmnts[middle].setExternalImpact
	                  (stratum->calculateImpactInPoint
								    (lmnts[middle].Cx, lmnts[middle].Cy));
	
	breaker.calculatePressure(&(lmnts[middle]), numOfCalcLmnts);
	
	calculateLmnts();
	if (fractionIsStopped) {
		numOfLmnts = numOfCalcLmnts;
		info("Fracture number", number, "is stopped.");
	}
	while ( (numOfCalcLmnts < numOfLmnts) ) {	
		double deltaBeta1 = calcAngleOfRotation(lmnts[front]);
		double deltaBeta2 = calcAngleOfRotation(lmnts[back]);
		
		// TODO - add the opportunity to work with a usual const tip element
		addNewLmnts(deltaBeta1, deltaBeta2);
		if (rotation == "predictor-corrector") {
			//	Clarifying the direction of fracture's growth like it is
			//	usually done in "predictor-corrector" methods
			calculateLmnts();
			deltaBeta1 = (calcAngleOfRotation(lmnts[front]) + deltaBeta1) / 2;
			deltaBeta2 = (calcAngleOfRotation(lmnts[back]) + deltaBeta2) / 2;		
			replaceTipElements(deltaBeta1, deltaBeta2);
		}
		calculateLmnts();
		if (fractionIsStopped) {
			numOfLmnts = numOfCalcLmnts;
			info("Fracture number", number, "is stopped.");
		}
	}
	numOfLmnts = numOfCalcLmnts;
	info("Fracture number", number, "is calculated.");
}

void Fracture::calculateLmnts() {
	double Ass, Asn, Ans, Ann;
	int N = 2 * numOfCalcLmnts;
	gsl_matrix *A = gsl_matrix_alloc(N, N);
	gsl_vector *b = gsl_vector_alloc(N);
	gsl_vector *x = gsl_vector_alloc(N);
	
	for (int i = 0; i < N; i += 2) {
		for (int j = 0; j < N; j += 2) {
			lmnts[front + j/2].calculateImpactOn(lmnts[front + i/2],
			                                     Ass, Asn, Ans, Ann);
			gsl_matrix_set(A, i, j, Ass);
			gsl_matrix_set(A, i, j + 1, Asn);
			gsl_matrix_set(A, i + 1, j, Ans);
			gsl_matrix_set(A, i + 1, j + 1, Ann);
		}
		gsl_vector_set(b, i, lmnts[front + i/2].getBs());
		gsl_vector_set(b, i + 1, lmnts[front + i/2].getBn());
	}

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
	
	for (int i = 0; i < N; i += 2) {
		lmnts[front + i/2].Ds = gsl_vector_get(x, i);
		lmnts[front + i/2].Dn = gsl_vector_get(x, i + 1);
	}
	
	gsl_matrix_free(A);
	gsl_vector_free(x);
	gsl_vector_free(b);
	
//	print(lmnts[front].K1() / sqrt(M_PI * halfLength) / sqrt(M_PI / 2), "\t",
//		lmnts[front].K2() / sqrt(M_PI * halfLength) / sqrt(M_PI / 2), "\t", 
//		lmnts[back].K1() / sqrt(M_PI * halfLength) / sqrt(M_PI / 2), "\t", 
//		lmnts[back].K2() / sqrt(M_PI * halfLength) / sqrt(M_PI / 2));
	
	if (lmnts[front].Dn > 0 || lmnts[back].Dn > 0)
		fractionIsStopped = true;
}

double Fracture::calcAngleOfRotation(const Element &lmnt1) const {
	double K1 = - lmnt1.Dn;
	double K2 = - lmnt1.Ds;
	Field tmp;
	double beta = 2 * tmp.arctan(- 2 * K2, K1 + sqrt(K1 * K1 + 8 * K2 * K2) );
	return beta;
}

void Fracture::addNewLmnts(const double &deltaBeta1, const double &deltaBeta2) {
	
	if (numOfCalcLmnts == 1) {
		addNewLmnts(deltaBeta1, deltaBeta2, 5 * a / 9);
		addNewLmnts(0, 0, a / 3);
		addNewLmnts(0, 0, a / 9);
	} else {
		double beta1 = lmnts[front].beta;
		double beta2 = lmnts[back].beta;
		front += 3;
		back -= 3;
		numOfCalcLmnts -= 6;
		halfLength -= 2 * a;
		double oldDeltaBeta1 = beta1 - lmnts[front].beta;
		double oldDeltaBeta2 = beta2 - lmnts[back].beta;
		
		addNewLmnts(oldDeltaBeta1, oldDeltaBeta2, a);
		addNewLmnts(deltaBeta1, deltaBeta2, 5 * a / 9);
		addNewLmnts(0, 0, a / 3);
		addNewLmnts(0, 0, a / 9);
	}
	
	breaker.calculatePressure(&lmnts[middle], numOfCalcLmnts);
}

void Fracture::addNewLmnts(const double& deltaBeta1, const double& deltaBeta2, 
                                                      const double& half_length) {
	double beta1 = lmnts[front].beta + deltaBeta1;
	double x1 = lmnts[front].Cx - lmnts[front].getA() * cos(lmnts[front].beta) 
	                            - half_length * cos(beta1);
	double y1 = lmnts[front].Cy - lmnts[front].getA() * sin(lmnts[front].beta)
	                            - half_length * sin(beta1);
	front--;
	lmnts[front] = Element(half_length, x1, y1, beta1, G, nu);
	lmnts[front].setExternalImpact(stratum->calculateImpactInPoint
	                                        (lmnts[front].Cx, lmnts[front].Cy));

	beta1 = lmnts[back].beta + deltaBeta2;
	x1 = lmnts[back].Cx + lmnts[back].getA() * cos(lmnts[back].beta) 
	                    + half_length * cos(beta1);
	y1 = lmnts[back].Cy + lmnts[back].getA() * sin(lmnts[back].beta)
	                    + half_length * sin(beta1);
	back++;
	lmnts[back] = Element(half_length, x1, y1, beta1, G, nu);
	lmnts[back].setExternalImpact(stratum->calculateImpactInPoint
	                                       (lmnts[back].Cx, lmnts[back].Cy));
	
	numOfCalcLmnts += 2;
	halfLength += 2 * half_length;
}

void Fracture::replaceTipElements(const double& deltaBeta1, const double& deltaBeta2) {
	back -= 3;
	front += 3;
	numOfCalcLmnts -= 6; 
	halfLength -= 2 * a;
	addNewLmnts(deltaBeta1, deltaBeta2, 5 * a / 9);
	addNewLmnts(0, 0, a / 3);
	addNewLmnts(0, 0, a / 9);
	breaker.calculatePressure(&lmnts[middle], numOfCalcLmnts);
}

int Fracture::getNumOfLmnts() const {
	return numOfCalcLmnts;
}

void Fracture::getPointsForPlot(double* x, double* y) const {
	x[0] = lmnts[front].Cx - lmnts[front].getA() * cos(lmnts[front].beta);
	y[0] = lmnts[front].Cy - lmnts[front].getA() * sin(lmnts[front].beta);
	
	for (int i = 0; i < numOfCalcLmnts; i++) {
		x[i+1] = lmnts[front + i].Cx + 
		         lmnts[front + i].getA() * cos(lmnts[front + i].beta);
		y[i+1] = lmnts[front + i].Cy + 
		         lmnts[front + i].getA() * sin(lmnts[front + i].beta);
	}
}

void Fracture::getPointsForDisplacementPlot(double* x, double* v) const {
	for (int i = 0; i < numOfCalcLmnts; i++) {
		x[i] = lmnts[front + i].Cx;
		v[i] = lmnts[front + i].Dn / 2;
	}
}

double Fracture::getHalfLength() const {
	return halfLength;
}
