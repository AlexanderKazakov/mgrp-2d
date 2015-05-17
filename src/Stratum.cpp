#include "Stratum.hpp"

Stratum::Stratum() {
}

Stratum::~Stratum() {
}

void Stratum::addFracture(int number, double volume, 
                          double x, double y, double beta, 
                          double halfLengthOfElements,
                          Breaker breaker) {
	info("Adding fracture number", number, "at (", x, ",", y, ") ...");
	fractures.push_back( Fracture(number, volume, halfLengthOfElements) );
	beginFracture = fractures.begin();
	fractures.back().allocateElements(x, y, beta, breaker);
}

void Stratum::setRheology(const double &_G, const double &_nu) {
	G = _G; nu = _nu;
	Element tmp;
	tmp.setRheology(G, nu);
}

void Stratum::setStresses(const double &_Sxx, const double &_Sxy, 
                          const double &_Syy) {
	Sxx = _Sxx;
	Sxy = _Sxy;
	Syy = _Syy;
}

void Stratum::setRanges(const double &_Xmin, const double &_Xmax,
                        const double &_Ymin, const double &_Ymax) {
	Xmin = _Xmin;
	Xmax = _Xmax;
	Ymin = _Ymin;
	Ymax = _Ymax;		
}

void Stratum::setSequence(const std::string _sequence) {
	sequence = _sequence;
}

void Stratum::setRotation(const std::string _rotation) {
	rotation = _rotation;
}

void Stratum::reserve(const int numberOfFractures) {
	fractures.reserve(numberOfFractures);
}

void Stratum::calculateTask() {
	info("Starting calculation of the task ...");
	
	if (sequence == "series") {
		beginFracture = fractures.begin();
		endFracture = fractures.begin() + 1;
		while (beginFracture != fractures.end()) {
			calculateStage();
			beginFracture->breakerIsInjecting = false;
			beginFracture++;
			if (endFracture != fractures.end()) endFracture++;
		}
		
	} else if (sequence == "parallel") {
		beginFracture = fractures.begin();
		endFracture = fractures.end();
		calculateStage();
		beginFracture = fractures.end();
		
	} else if (sequence == "series with feedback") {
		beginFracture = fractures.begin();
		endFracture = fractures.begin();
		do {
			endFracture++;
			calculateStage();
			(endFracture - 1)->breakerIsInjecting = false;
		} while (endFracture != fractures.end());
		beginFracture = fractures.end();
	}
	
	info("Calculation of the task is done.");
}

Field Stratum::calculateImpactInPoint(const double& x, const double& y) const{
	auto fracture = fractures.begin();
	Field field;
	while (fracture != beginFracture) {
		field += fracture->calculateImpactInPoint(x, y);
		fracture++;
	}
	field.Sxx += Sxx;
	field.Sxy += Sxy;
	field.Syy += Syy;
	return field;
}

void Stratum::calculateStage() {
	info("Starting calculation of the next stage. Fractures from", 
	     beginFracture->getNumber(), "to", 
	     (endFracture - 1)->getNumber(), "...");
	
	for(auto frac1 = beginFracture; frac1 != endFracture; frac1++)
		frac1->setExternalImpactAndBreakerPressure();
	calculateElements();
	bool calculationIsCompleted = true;
	for(auto frac1 = beginFracture; frac1 != endFracture; frac1++)
		calculationIsCompleted *= frac1->isCompleted();

	while( ! calculationIsCompleted ) {
		for(auto frac1 = beginFracture; frac1 != endFracture; frac1++)
			frac1->grow();
		if (rotation == "predictor-corrector") {
			// Clarifying the direction of fracture's growth like it is
			// usually done in "predictor-corrector" methods
			calculateElements();
			for(auto frac1 = beginFracture; frac1 != endFracture; frac1++)
				frac1->correctRotation();
		}
		calculateElements();
		calculationIsCompleted = true;
		for(auto frac1 = beginFracture; frac1 != endFracture; frac1++)
			calculationIsCompleted *= frac1->isCompleted();
	}
	info("Calculation of stage is complete.");
}

void Stratum::calculateElements() {
	int N = 0;
	for(auto frac1 = beginFracture; frac1 != endFracture; frac1++) {
		N += 2 * frac1->getNumOfElements();
	}
	gsl_matrix *A = gsl_matrix_alloc(N, N);
	gsl_vector *b = gsl_vector_alloc(N);
	gsl_vector *x = gsl_vector_alloc(N);

	int i1 = 0;
	for(auto frac1 = beginFracture; frac1 != endFracture; frac1++) {
		int i2 = 0;
		for(auto frac2 = beginFracture; frac2 != endFracture; frac2++) {
			frac1->fillInMatrixA(*frac2, A, i1, i2);
			i2 += 2 * frac2->getNumOfElements();
		}
		frac1->fillInVectorB(b, i1);
		i1 += 2 * frac1->getNumOfElements();
	}
	
	// to save original of A after LU-decomposition
	gsl_matrix *LU = gsl_matrix_alloc(N, N);
	gsl_matrix_memcpy(LU, A);	
	gsl_permutation *p = gsl_permutation_alloc(N);
	int signum;		
	gsl_linalg_LU_decomp(LU, p, &signum);
	gsl_linalg_LU_solve(LU, p, b, x);
	gsl_permutation_free(p);
	gsl_matrix_free(LU);
	
	checkSLE(A, x, b, N);
		
	int i = 0;
	for(auto frac1 = beginFracture; frac1 != endFracture; frac1++) {
		frac1->takeDDfromVectorX(x, i);
		i += 2 * frac1->getNumOfElements();
	}
	
	gsl_matrix_free(A);
	gsl_vector_free(x);
	gsl_vector_free(b);
}


void Stratum::visualize() const{
	info("Starting visualisation ...");
	mglGraph gr = mglGraph(0, 1200, 800);
	gr.SetRanges(Xmin, Xmax, Ymin, Ymax);
	gr.Axis();

	drawFractures(gr);
//	drawField(gr);
//	drawStressDirections(gr);
	
	gr.WriteFrame("fractures.png");
	info("Visualisation is done.");
}

void Stratum::drawDisplacements() const {
	info("Drawing displacements ...");
	auto fracture = fractures.begin();
	int N = fracture->getNumOfElements();
	double *_x = new double[N];
	double *_v = new double[N];
	// numerical results
	fracture->getPointsForDisplacementPlot(_x, _v);
	mglData v;
	mglData Cx;
	Cx.Set(_x, N);
	v.Set(_v, N);
	
	mglGraph gr = mglGraph(0, 1200, 800);
	gr.SetRanges(1.1*_x[0], 1.1*_x[N-1], 1.1*_v[N/2 + 1], - 1.1*_v[N/2 + 1]);
	gr.Axis();

	gr.Plot(Cx, v, ".k");
	for (int i = 0; i < N; i++) {
		v.a[i] = -v.a[i];
	}
	gr.Plot(Cx, v, ".k");
	
	std::string pressureType;
	double _a, b, c;
	fractures.begin()->getBreaker()->getType(_a, b, c, pressureType);
	// analytical solution
	if ( pressureType == "const" ) {
		for (int i = 0; i < N; i++) {
			double l = fracture->getLeftLength();
			double x = Cx.a[i];
			v.a[i] = (Syy - c) / G * (1 - nu) * sqrt(l*l - x*x);
		}
	} else if ( pressureType == "lag" ) {
		for (int i = 0; i < N; i++) {
			double l = fracture->getLeftLength();
			double x = Cx.a[i];
			double p = - c;
			double sigma = Syy - c;
			double a = (1 - _a) * l;
			v.a[i] =  sigma / G * (1 - nu) * sqrt(l*l - x*x)
					
			        - 2 * p * (1 - nu) / G / M_PI *
			            ( acos((l - a) / l) * sqrt(l*l - x*x) 
					
			            + x / 2 * log(fabs(
			                (x * sqrt(2*a*l - a*a) + (l - a) * sqrt(l*l - x*x)) /
			                ( x * sqrt(2*a*l - a*a) - (l - a) *	sqrt(l*l - x*x)) ))
			        
			            - (l - a) / 2 * log(fabs(  
			                (sqrt(2*a*l - a*a) + sqrt(l*l - x*x)) / 
			                (sqrt(2*a*l - a*a) - sqrt(l*l - x*x)) )) );
		}
	}
	gr.Plot(Cx, v, "r");
	for (int i = 0; i < N; i++) {
		v.a[i] = -v.a[i];
	}
	gr.Plot(Cx, v, "r");
	gr.WriteFrame("displacements.png");
	delete [] _x;
	delete [] _v;
}

void Stratum::drawFractures(mglGraph& gr) const{
	info("Drawing fractures ...");
	auto fracture = fractures.begin();
	while (fracture != fractures.end()) {
		int N = fracture->getNumOfElements() + 1;
		double *_x = new double[N];
		double *_y = new double[N];
		fracture->getPointsForPlot(_x, _y);
		mglData x;
		mglData y;
		x.Set(_x, N);
		y.Set(_y, N);
		gr.Plot(x, y, ".g");
		delete [] _x;
		delete [] _y;
		fracture++;
	}
}

void Stratum::drawField(mglGraph &gr) const{
	info("Drawing field ...");
	int N = 104;
	mglData x(N);
	mglData y(N);
	mglData f(N, N);
	double maxF = 0;
	for (int i = 0; i < N; i += 1) {
		x.a[i] = (Xmax - Xmin) * i / (N-1) + Xmin;
		y.a[i] = (Ymax - Ymin) * i / (N-1) + Ymin;
	}
	for (int i = 0; i < N; i += 1)
		for (int j = 0; j < N; j += 1) {
			double _x = x.a[i];
			double _y = y.a[j];
			double _f = calculateImpactInPoint(_x, _y).Trace();
			maxF = (maxF > fabs(_f)) ? maxF : fabs(_f);
			f.a[i + N * j] = _f;
		}
	maxF = maxF / 10;
	for (int i = 0; i < N; i += 1)
		for (int j = 0; j < N; j += 1) {
			f.a[i + N * j] = f.a[i + N * j] / maxF;
		}
	gr.Dens(x, y, f);
	gr.Colorbar();
}

void Stratum::drawStressDirections(mglGraph &gr) const{
	info("Drawing main stress directions ...");
	int N = 23;
	mglData x(N);
	mglData y(N);
	mglData ax(N, N);
	mglData ay(N, N);
	for (int i = 1; i < N; i += 1) {
		x.a[i] = (Xmax - Xmin) * i / N + Xmin;
		y.a[i] = (Ymax - Ymin) * i / N + Ymin;
	}
	for (int i = 1; i < N; i += 1)
		for (int j = 1; j < N; j += 1) {
			double _x = x.a[i];
			double _y = y.a[j];
			Field localField = calculateImpactInPoint(_x, _y);
			double Smax = localField.Smax() * 5;			
			ax.a[i + N * j] = Smax * cos(localField.directionOfMaxTensileStress());
			ay.a[i + N * j] = Smax * sin(localField.directionOfMaxTensileStress());
		}
	gr.Vect(x, y, ax, ay, "0");
	for (int i = 1; i < N; i += 1)
		for (int j = 1; j < N; j += 1) {
			ax.a[i + N * j] = - ax.a[i + N * j];
			ay.a[i + N * j] = - ay.a[i + N * j];
		}
	gr.Vect(x, y, ax, ay, "0");
}
