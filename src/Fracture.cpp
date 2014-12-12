#include "Fracture.hpp"

Fracture::Fracture() {
	breaks = NULL;
}

Fracture::Fracture(Stratum *stratum, int number, double h_length, int numOfElms,
		double a, double b, double c, std::string pressureType, std::string tip,
		std::string rotation): stratum(stratum), number(number), tip(tip),
		half_lengthOfBreaks(h_length), numOfBrks(numOfElms), rotation(rotation) {
	breaks = NULL;
	fractionIsStopped = false;
	double _G, _nu;
	stratum->getRheology(_G, _nu);
	G = _G;	nu = _nu;
	front = back = middle = numOfBrks/2;
	fluid.setType(a, b, c, pressureType);
	numOfCalcBrks = 1;
}

Fracture::~Fracture() {
	if (breaks) {
		// Deletion of dynamical allocated memory (!)
		delete [] breaks;
	}
}

void Fracture::allocateBreaks(double x, double y, double beta) {
	// Dynamical memory allocation (!)
	breaks = new Break[numOfBrks];
	breaks[middle] = Break(half_lengthOfBreaks, x, y, beta, G, nu);
}

Field Fracture::calculateImpactInPoint(const double &x, const double &y) const {
	Field field;
	for (int i = front; i <= back; i++)
		field += breaks[i].calculateImpactInPoint(x, y);
	return field;
}

void Fracture::calculate() {
	breaks[middle].setExternalImpact(stratum->calculateImpactInPoint
								(breaks[middle].Cx, breaks[middle].Cy));

	fluid.calculatePressure(&(breaks[middle]), numOfCalcBrks);
	
	calculateBreaks();
	if (fractionIsStopped) {
		numOfBrks = numOfCalcBrks;
		std::cout << "Fracture number " << number << " is stopped!\n";
	}
	while ( (numOfCalcBrks < numOfBrks) ) {	
		double deltaBeta1 = calcAngleOfRotation(breaks[front]);
		double deltaBeta2 = calcAngleOfRotation(breaks[back]);
		
		addNewBreaks(deltaBeta1, deltaBeta2);
		if (rotation == "predictor-corrector") {
			//	Clarifying the direction of fracture's growth like it is
			//	usually done in "predictor-corrector" methods
			calculateBreaks();
			deltaBeta1 = (calcAngleOfRotation(breaks[front]) + deltaBeta1) / 2;
			deltaBeta2 = (calcAngleOfRotation(breaks[back]) + deltaBeta2) / 2;
		
			back--; front++; numOfCalcBrks -= 2;
			addNewBreaks(deltaBeta1, deltaBeta2);
		}
		calculateBreaks();
		if (fractionIsStopped) {
			numOfBrks = numOfCalcBrks;
			std::cout << "Fracture number " << number << " is stopped!\n";
		}
	}
}

void Fracture::calculateBreaks() {
	double Ass, Asn, Ans, Ann;
	int N = 2 * numOfCalcBrks;
	gsl_matrix *A = gsl_matrix_alloc(N, N);
	gsl_vector *b = gsl_vector_alloc(N);
	gsl_vector *x = gsl_vector_alloc(N);
	
	for (int i = 0; i < N; i += 2) {
		for (int j = 0; j < N; j += 2) {
			breaks[front + j/2].calculateImpactOn(breaks[front + i/2],
													Ass, Asn, Ans, Ann);
			gsl_matrix_set(A, i, j, Ass);
			gsl_matrix_set(A, i, j + 1, Asn);
			gsl_matrix_set(A, i + 1, j, Ans);
			gsl_matrix_set(A, i + 1, j + 1, Ann);
		}
		gsl_vector_set(b, i, breaks[front + i/2].getBs());
		gsl_vector_set(b, i + 1, breaks[front + i/2].getBn());
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
		breaks[front + i/2].Ds = gsl_vector_get(x, i);
		breaks[front + i/2].Dn = gsl_vector_get(x, i + 1);
	}
	
	gsl_matrix_free(A);
	gsl_vector_free(x);
	gsl_vector_free(b);
	
	if (breaks[front].Dn > 0 || breaks[back].Dn > 0)
		fractionIsStopped = true;
}

double Fracture::calcAngleOfRotation(const Break &break1) const {
	double K1 = - break1.Dn;
	double K2 = - break1.Ds;
	Field tmp;
	double beta = 2 * tmp.arctan(- 2 * K2, K1 + sqrt(K1 * K1 + 8 * K2 * K2) );
	return beta;
}

void Fracture::addNewBreaks(const double &deltaBeta1, const double &deltaBeta2) {
	double beta1 = breaks[front].beta + deltaBeta1;
	double x1 = breaks[front].Cx - half_lengthOfBreaks *
			(cos(beta1) + cos(breaks[front].beta));
	double y1 = breaks[front].Cy - half_lengthOfBreaks *
			(sin(beta1) + sin(breaks[front].beta));
	front--;
	breaks[front] = Break(half_lengthOfBreaks, x1, y1, beta1, G, nu);
	breaks[front].setExternalImpact(stratum->calculateImpactInPoint
			(breaks[front].Cx, breaks[front].Cy));

	beta1 = breaks[back].beta + deltaBeta2;
	x1 = breaks[back].Cx + half_lengthOfBreaks *
			(cos(beta1) + cos(breaks[back].beta));
	y1 = breaks[back].Cy + half_lengthOfBreaks *
			(sin(beta1) + sin(breaks[back].beta));
	back++;
	breaks[back] = Break(half_lengthOfBreaks, x1, y1, beta1, G, nu);
	breaks[back].setExternalImpact(stratum->calculateImpactInPoint
			(breaks[back].Cx, breaks[back].Cy));
	
	numOfCalcBrks += 2;
	fluid.calculatePressure(&breaks[middle], numOfCalcBrks);
}

int Fracture::getNumOfBreaks() const {
	return numOfBrks;
}

void Fracture::getPointsForPlot(double* x, double* y) const {
	x[0] = breaks[front].Cx - half_lengthOfBreaks * cos(breaks[front].beta);
	y[0] = breaks[front].Cy - half_lengthOfBreaks * sin(breaks[front].beta);
	
	for (int i = 0; i < numOfCalcBrks; i++) {
		x[i+1] = breaks[front + i].Cx + 
				half_lengthOfBreaks * cos(breaks[front + i].beta);
		y[i+1] = breaks[front + i].Cy + 
				half_lengthOfBreaks * sin(breaks[front + i].beta);
	}
}
