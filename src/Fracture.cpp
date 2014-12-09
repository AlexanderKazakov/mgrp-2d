#include "Fracture.hpp"

Fracture::Fracture() {
	breaks = NULL;
}

Fracture::Fracture(Stratum *stratum, int number, double h_length, int numOfElms,
		double a, double b, double c, std::string pressureType): stratum(stratum),
		number(number), half_lengthOfBreaks(h_length), numOfBrks(numOfElms) {
	breaks = NULL;
	double _G, _nu;
	stratum->getRheology(_G, _nu);
	G = _G;	nu = _nu;
	front = back = middle = numOfBrks/2;
	fluid.setType(a, b, c, pressureType);
	numOfCalcBrks = 1;
}

Fracture::~Fracture() {
	if (breaks) {
		delete [] breaks;
	}
}

void Fracture::allocateBreaks(double x, double y, double beta) {
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
								(breaks[middle].getCx(), breaks[middle].getCy()));

	fluid.calculatePressure(&(breaks[middle]), numOfCalcBrks);

	//	While the fracture is not stopped
	while ( (numOfCalcBrks < numOfBrks) ) {
		//	Calculate the existing configuration
		bool fractionIsStopped = calculateBreaks();
		if (fractionIsStopped) {
			numOfBrks = numOfCalcBrks;
			std::cout << "Fracture number " << number << " is stopped!\n";
		}
		else {
			Break *break1 = &(breaks[front]);
			//	Determine the direction of the fraction's growth
			double K1 = - break1->getDn();
			double K2 = - break1->getDs();
			double beta1 = break1->getBeta() + calcAngleOfRotation(K1, K2);
			double x1 = break1->getCx() - half_lengthOfBreaks * 
						( cos(beta1) + cos(break1->getBeta()) );
			double y1 = break1->getCy() - half_lengthOfBreaks *
						( sin(beta1) + sin(break1->getBeta()) );
			front--;
			break1 = &(breaks[front]);
			*break1 = Break(half_lengthOfBreaks, x1, y1, beta1, G, nu);
			break1->setExternalImpact(stratum->calculateImpactInPoint
											(break1->getCx(), break1->getCy()));

			break1 = &(breaks[back]);
			K1 = - break1->getDn();
			K2 = - break1->getDs();
			beta1 = break1->getBeta() + calcAngleOfRotation(K1, K2);
			x1 = break1->getCx() + half_lengthOfBreaks *
					(cos(beta1) + cos(break1->getBeta()));
			y1 = break1->getCy() + half_lengthOfBreaks *
					(sin(beta1) + sin(break1->getBeta()));
			back++;
			break1 = &(breaks[back]);
			*break1 = Break(half_lengthOfBreaks, x1, y1, beta1, G, nu);
			break1->setExternalImpact(stratum->calculateImpactInPoint
					(break1->getCx(), break1->getCy()));

			numOfCalcBrks += 2;
			fluid.calculatePressure(&breaks[middle], numOfCalcBrks);
		}
//		for (int i = 0; i < numOfBrks; i++) {
//			std::cout << breaks[i] << std::endl;
//		}
//		std::cout << std::endl;
	}
	calculateBreaks();
}

bool Fracture::calculateBreaks() {
	double Ass, Asn, Ans, Ann;
	int N = 2 * numOfCalcBrks;
	gsl_matrix *A = gsl_matrix_alloc(N, N);
	gsl_vector *b = gsl_vector_alloc(N);
	gsl_vector *x = gsl_vector_alloc(N);
	
	for (int i = 0; i < N; i += 2) {
		Break *break1 = &(breaks[front + i/2]);
		for (int j = 0; j < N; j += 2) {
			Break *break2 = &(breaks[front + j/2]);
			break1->calculateImpactOf(*break2, Ass, Asn, Ans, Ann);
			gsl_matrix_set(A, i, j, Ass);
			gsl_matrix_set(A, i, j + 1, Asn);
			gsl_matrix_set(A, i + 1, j, Ans);
			gsl_matrix_set(A, i + 1, j + 1, Ann);
		}
		gsl_vector_set(b, i, break1->getBs());
		gsl_vector_set(b, i + 1, break1->getBn());
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
		Break *break1 = &(breaks[front + i/2]);
		break1->setDs(gsl_vector_get(x, i));
		break1->setDn(gsl_vector_get(x, i + 1));
	}
	
	gsl_matrix_free(A);
	gsl_vector_free(x);
	gsl_vector_free(b);
	
	bool fractionIsStopped = false;
	if (breaks[front].getDn() > 0 || 
				breaks[back].getDn() > 0)
		fractionIsStopped = true;
	
	return fractionIsStopped;
}

double Fracture::calcAngleOfRotation(const double& K1, const double& K2) const {
	Field tmp;
	double beta = 2 * tmp.arctan(- 2 * K2, K1 + sqrt(K1 * K1 + 8 * K2 * K2) );
	return beta;
}

int Fracture::getNumOfBreaks() const {
	return numOfBrks;
}

void Fracture::getPointsForPlot(double* x, double* y) const {
	Break *break1 = &(breaks[front]);
		x[0] = break1->getCx() - half_lengthOfBreaks * cos(break1->getBeta());
		y[0] = break1->getCy() - half_lengthOfBreaks * sin(break1->getBeta());
	
	for (int i = 0; i < numOfCalcBrks; i++) {
		Break *break1 = &(breaks[front + i]);
		x[i+1] = break1->getCx() + half_lengthOfBreaks * cos(break1->getBeta());
		y[i+1] = break1->getCy() + half_lengthOfBreaks * sin(break1->getBeta());
	}
}

bool Fracture::operator==(const Fracture &other) const {
	return ( number == other.getNumber() );
}

bool Fracture::operator!=(const Fracture &other) const {
	return ( number != other.getNumber() );
}	

bool Fracture::operator<(const Fracture &other) const {
	return ( number < other.getNumber() );
}

int Fracture::getNumber() const {
	return number;
}