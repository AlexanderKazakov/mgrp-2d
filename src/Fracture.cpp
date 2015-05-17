#include "Fracture.hpp"

Stratum *Fracture::stratum;
std::string Fracture::tip;

Fracture::Fracture() {
	breaker = NULL;
}

Fracture::Fracture(int number, double volume, double halfLengthOfElements):
                   number(number), taskVolume(volume), a(halfLengthOfElements),
                   fractionIsStoppedL(false),
                   fractionIsStoppedR(false),
                   breakerIsInjecting(true) {
	breaker = NULL;
}

Fracture::~Fracture() {
	if (breaker)
		delete breaker;
}

void Fracture::allocateElements(double x, double y, double beta, Breaker _breaker) {
	breaker = new Breaker;
	*breaker = _breaker; 
	
	elementsL.reserve((int) (5 * volume / a));
	elementsR.reserve((int) (5 * volume / a));
	
	if (tip == "135") {
		// Add three small tip elements to the corners
		elementsR.push_back(Element(5 * a / 9, x + 5 * a / 9 * cos(beta), 
											   y + 5 * a / 9 * sin(beta), beta));
		elementsR.push_back(Element(a / 3, x + 13 * a / 9 * cos(beta), 
										   y + 13 * a / 9 * sin(beta), beta));
		elementsR.push_back(Element(a / 9, x + 17 * a / 9 * cos(beta), 
										   y + 17 * a / 9 * sin(beta), beta));
		elementsL.push_back(Element(5 * a / 9, x - 5 * a / 9 * cos(beta), 
											   y - 5 * a / 9 * sin(beta), beta));
		elementsL.push_back(Element(a / 3, x - 13 * a / 9 * cos(beta), 
										   y - 13 * a / 9 * sin(beta), beta));
		elementsL.push_back(Element(a / 9, x - 17 * a / 9 * cos(beta), 
										   y - 17 * a / 9 * sin(beta), beta));
	} else if (tip == "const") {
		// Add one usual element to the corners
		elementsL.push_back(Element(a, x - a * cos(beta), y - a * sin(beta), beta));
		elementsR.push_back(Element(a, x + a * cos(beta), y + a * sin(beta), beta));
	}
	
	leftLength = rightLength = 2*a;
	
	// Uncomment for creating a special initial shapes of the fracture 
	createLinearShape(1.0);
	grow(M_PI/36, M_PI/36);
}

bool Fracture::isCompleted() const {
	return (volume > taskVolume) || isStopped() || ( ! breakerIsInjecting );
}

bool Fracture::isStopped() const {
	return (fractionIsStoppedL && fractionIsStoppedR);
}

Field Fracture::calculateImpactInPoint(const double &x, const double &y) const {
	Field field;
	for (auto element = elementsL.begin(); element != elementsL.end(); element++) {
		field += element->calculateImpactInPoint(x, y);
	}
	for (auto element = elementsR.begin(); element != elementsR.end(); element++) {
		field += element->calculateImpactInPoint(x, y);
	}
	return field;
}

void Fracture::grow() {
	double deltaBeta1 = elementsL.back().calcAngleOfRotation();
	double deltaBeta2 = elementsR.back().calcAngleOfRotation();
	print("predictor ", "K1 = ", elementsR.back().K1(), "\tK2 = ", elementsR.back().K2());
	deltaBetaL = deltaBeta1; deltaBetaR = deltaBeta2;
	// Fracture grows after injecting
	//if ( ! isStopped() ) grow(deltaBeta1, deltaBeta2);
	// Fracture doesn't grow after injecting
	if ( ! isCompleted() ) grow(deltaBeta1, deltaBeta2);
}

void Fracture::correctRotation() {
	double deltaBeta1 = (elementsL.back().calcAngleOfRotation() + deltaBetaL) / 2;
	double deltaBeta2 = (elementsR.back().calcAngleOfRotation() + deltaBetaR) / 2;
	print("corrector ", "K1 = ", elementsR.back().K1(), "\tK2 = ", elementsR.back().K2());
	if ( ! isStopped() ) replaceTipElements(deltaBeta1, deltaBeta2);
}

void Fracture::setExternalImpactAndBreakerPressure() {
	for (auto element = elementsL.begin(); element != elementsL.end(); element++) 
		element->setExternalImpact
	                     (stratum->calculateImpactInPoint(element->Cx, element->Cy));
	for (auto element = elementsR.begin(); element != elementsR.end(); element++) 
		element->setExternalImpact
	                     (stratum->calculateImpactInPoint(element->Cx, element->Cy));
	breaker->calculatePressure(this);
}


void Fracture::fillInMatrixA(const Fracture& frac2, gsl_matrix* A, 
                                        int i1, int i2) const {
	double Ass, Asn, Ans, Ann;
	int i = i1;
	for (auto element1 = elementsL.rbegin(); element1 != elementsL.rend(); element1++) {
		int j = i2;
		for (auto element2 = frac2.elementsL.rbegin(); 
		          element2 != frac2.elementsL.rend(); element2++) {
			element2->calculateImpactOn(*element1, Ass, Asn, Ans, Ann);
			gsl_matrix_set(A, i, j, Ass);
			gsl_matrix_set(A, i, j + 1, Asn);
			gsl_matrix_set(A, i + 1, j, Ans);
			gsl_matrix_set(A, i + 1, j + 1, Ann);
			j += 2;
		}
		for (auto element2 = frac2.elementsR.begin(); 
		          element2 != frac2.elementsR.end(); element2++) {
			element2->calculateImpactOn(*element1, Ass, Asn, Ans, Ann);
			gsl_matrix_set(A, i, j, Asn);
			gsl_matrix_set(A, i, j + 1, Ass);
			gsl_matrix_set(A, i + 1, j, Ann);
			gsl_matrix_set(A, i + 1, j + 1, Ans);
			j += 2;
		}
		i += 2;
	}
	for (auto element1 = elementsR.begin(); element1 != elementsR.end(); element1++) {
		int j = i2;
		for (auto element2 = frac2.elementsL.rbegin(); 
		          element2 != frac2.elementsL.rend(); element2++) {
			element2->calculateImpactOn(*element1, Ass, Asn, Ans, Ann);
			gsl_matrix_set(A, i, j, Ans);
			gsl_matrix_set(A, i, j + 1, Ann);
			gsl_matrix_set(A, i + 1, j, Ass);
			gsl_matrix_set(A, i + 1, j + 1, Asn);
			j += 2;
		}
		for (auto element2 = frac2.elementsR.begin(); 
		          element2 != frac2.elementsR.end(); element2++) {
			element2->calculateImpactOn(*element1, Ass, Asn, Ans, Ann);
			gsl_matrix_set(A, i, j, Ann);
			gsl_matrix_set(A, i, j + 1, Ans);
			gsl_matrix_set(A, i + 1, j, Asn);
			gsl_matrix_set(A, i + 1, j + 1, Ass);
			j += 2;
		}
		i += 2;
	}
}

void Fracture::fillInVectorB(gsl_vector* b, int i) const {
	for (auto element1 = elementsL.rbegin(); element1 != elementsL.rend(); element1++) {
		gsl_vector_set(b, i, element1->getBs());
		gsl_vector_set(b, i + 1, element1->getBn());
		i += 2;
	}
	for (auto element1 = elementsR.begin(); element1 != elementsR.end(); element1++) {
		gsl_vector_set(b, i + 1, element1->getBs());
		gsl_vector_set(b, i, element1->getBn());
		i += 2;
	}
}

void Fracture::takeDDfromVectorX(const gsl_vector* x, int i) {
	volume = 0;
	for (auto element1 = elementsL.rbegin(); element1 != elementsL.rend(); element1++) {
		element1->Ds = gsl_vector_get(x, i);
		element1->Dn = gsl_vector_get(x, i + 1);
		if (element1->Dn > 0)
			element1->Dn = 0;
		volume += - element1->Dn * 2 * element1->a;
		i += 2;
	}
	for (auto element1 = elementsR.begin(); element1 != elementsR.end(); element1++) {
		element1->Ds = gsl_vector_get(x, i + 1);
		element1->Dn = gsl_vector_get(x, i);
		if (element1->Dn > 0)
			element1->Dn = 0;
		volume += - element1->Dn * 2 * element1->a;
		i += 2;
	}
	if (elementsL.back().Dn == 0) {
		fractionIsStoppedL = true;
		info("Fracture number", number, "is stopped at the left corner");
	}
	if (elementsR.back().Dn == 0) {
		fractionIsStoppedR = true;
		info("Fracture number", number, "is stopped at the right corner");
	}
}

int Fracture::getNumber() const {
	return number;
}

int Fracture::getNumOfElements() const {
	return elementsL.size() + elementsR.size();
}

void Fracture::getPointsForPlot(double* x, double* y) const {
	x[0] = elementsL.back().Cx - elementsL.back().a * cos(elementsL.back().beta);
	y[0] = elementsL.back().Cy - elementsL.back().a * sin(elementsL.back().beta);

	int i = 1;
	for (auto element = elementsL.rbegin(); element != elementsL.rend(); element++) {
		x[i] = element->Cx + element->a * cos(element->beta);
		y[i] = element->Cy + element->a * sin(element->beta);
		i++;
	}
	for (auto element = elementsR.begin(); element != elementsR.end(); element++) {
		x[i] = element->Cx + element->a * cos(element->beta);
		y[i] = element->Cy + element->a * sin(element->beta);
		i++;
	}
}

void Fracture::getPointsForDisplacementPlot(double* x, double* v) const {
	int i = 0;
	for (auto element = elementsL.rbegin(); element != elementsL.rend(); element++) {
		x[i] = element->Cx;
		v[i] = element->Dn / 2;
		i++;
	}
	for (auto element = elementsR.begin(); element != elementsR.end(); element++) {
		x[i] = element->Cx;
		v[i] = element->Dn / 2;
		i++;
	}
}

double Fracture::getLeftLength() const {
	return leftLength;
}

double Fracture::getRightLength() const {
	return rightLength;
}

Breaker *Fracture::getBreaker() const {
	return breaker;
}

void Fracture::setStratum(Stratum* _stratum){
	stratum = _stratum;
}

void Fracture::setTip(std::string _tip){
	tip = _tip;
}

void Fracture::grow(const double &deltaBeta1, const double &deltaBeta2) {
	if (tip == "135") {
		if ( (elementsL.size() == 3) && (elementsR.size() == 3) ) {
			double beta0 = elementsL.back().beta;
			double x0 = elementsL.front().Cx + 5 * a / 9 * cos(beta0);
			double y0 = elementsL.front().Cy + 5 * a / 9 * sin(beta0);
			if ( ! fractionIsStoppedL ) {
				elementsL.pop_back(); elementsL.pop_back(); elementsL.pop_back();
				elementsL.push_back(Element(a, x0 - a * cos(beta0), 
											   y0 - a * sin(beta0), beta0));
				elementsL.back().setExternalImpact
									 (stratum->calculateImpactInPoint
										  (elementsL.back().Cx, elementsL.back().Cy));
			}
			if ( ! fractionIsStoppedR ) {
				elementsR.pop_back(); elementsR.pop_back(); elementsR.pop_back();
				elementsR.push_back(Element(a, x0 + a * cos(beta0), 
											   y0 + a * sin(beta0), beta0));
				elementsR.back().setExternalImpact
									 (stratum->calculateImpactInPoint
										  (elementsR.back().Cx, elementsR.back().Cy));
			}
			addNewElements(deltaBeta1, deltaBeta2, 5 * a / 9);
			addNewElements(0, 0, a / 3);
			addNewElements(0, 0, a / 9);
		} else {
			double beta1 = elementsL.back().beta;
			double beta2 = elementsR.back().beta;
			if ( ! fractionIsStoppedL ) {
				elementsL.pop_back(); elementsL.pop_back(); elementsL.pop_back();
				leftLength -= 2*a;
			}
			if ( ! fractionIsStoppedR ) {
				elementsR.pop_back(); elementsR.pop_back(); elementsR.pop_back();
				rightLength -= 2*a;
			}
			double oldDeltaBeta1 = beta1 - elementsL.back().beta;
			double oldDeltaBeta2 = beta2 - elementsR.back().beta;

			addNewElements(oldDeltaBeta1, oldDeltaBeta2, a);
			addNewElements(deltaBeta1, deltaBeta2, 5 * a / 9);
			addNewElements(0, 0, a / 3);
			addNewElements(0, 0, a / 9);
		}
	} else if (tip == "const") {
		addNewElements(deltaBeta1, deltaBeta2, a);
	}
	breaker->calculatePressure(this);
}

void Fracture::addNewElements(const double& deltaBeta1, const double& deltaBeta2, 
                              const double& half_length) {
	if ( ! fractionIsStoppedL ) {
		double beta1 = elementsL.back().beta + deltaBeta1;
		double x1 = elementsL.back().Cx 
		          - elementsL.back().a * cos(elementsL.back().beta) 
		          - half_length * cos(beta1);
		double y1 = elementsL.back().Cy 
		          - elementsL.back().a * sin(elementsL.back().beta)
		          - half_length * sin(beta1);

		elementsL.push_back(Element(half_length, x1, y1, beta1));
		elementsL.back().setExternalImpact
		                     (stratum->calculateImpactInPoint
		                          (elementsL.back().Cx, elementsL.back().Cy));
		leftLength += 2 * half_length;
	}

	if ( ! fractionIsStoppedR ) {
		double beta1 = elementsR.back().beta + deltaBeta2;
		double x1 = elementsR.back().Cx 
		          + elementsR.back().a * cos(elementsR.back().beta) 
		          + half_length * cos(beta1);
		double y1 = elementsR.back().Cy 
		          + elementsR.back().a * sin(elementsR.back().beta)
		          + half_length * sin(beta1);

		elementsR.push_back(Element(half_length, x1, y1, beta1));
		elementsR.back().setExternalImpact
		                     (stratum->calculateImpactInPoint
		                           (elementsR.back().Cx, elementsR.back().Cy));
		rightLength += 2 * half_length;
	}
}

void Fracture::replaceTipElements(const double& deltaBeta1, const double& deltaBeta2) {
	if (tip == "135") {
		if ( ! fractionIsStoppedL ) {
			elementsL.pop_back(); elementsL.pop_back(); elementsL.pop_back();
			leftLength -= 2*a;
		}
		if ( ! fractionIsStoppedR ) {
			elementsR.pop_back(); elementsR.pop_back(); elementsR.pop_back();
			rightLength -= 2*a;
		}
		addNewElements(deltaBeta1, deltaBeta2, 5 * a / 9);
		addNewElements(0, 0, a / 3);
		addNewElements(0, 0, a / 9);
	} else if (tip == "const") {
		if ( ! fractionIsStoppedL ) {
			elementsL.pop_back();
			leftLength -= 2*a;
		}
		if ( ! fractionIsStoppedR ) {
			elementsR.pop_back();
			rightLength -= 2*a;
		}
		addNewElements(deltaBeta1, deltaBeta2, a);
	}
	breaker->calculatePressure(this);
}

void Fracture::createLinearShape(const double& initHalfLength) {
	while (leftLength < initHalfLength - 2*a) {
		grow(0.0, 0.0);
	}
}
