#include "Stratum.hpp"

Stratum::Stratum() {
	info("Stratum is created.");
}

Stratum::~Stratum() {
}

void Stratum::addFracture(int number, int numOfLmnts, 
                          double x, double y, double beta, 
                          double halfLengthOfLmnts,
                          double a, double b, double c,
                          std::string pressureType, std::string tip, 
                          std::string rotation) {
	info("Adding fracture number", number, "at (", x, ",", y, ") ...");
	fractures.push_back( Fracture(this, number, numOfLmnts, halfLengthOfLmnts, 
                                  a, b, c, pressureType, tip, rotation));
	fractures.back().allocateLmnts(x, y, beta);
	breakerOfFirstFracture.setType(a, b, c, pressureType);
}

void Stratum::setRheology(const double &_G, const double &_nu) {
	G = _G; nu = _nu;
}

void Stratum::getRheology(double &_G, double &_nu) const {
	_G = G; _nu = nu;
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

void Stratum::reserve(const int numberOfFractures) {
	fractures.reserve(numberOfFractures);
}

void Stratum::calculate() {
	info("Starting calculation ...");
	currentFracture = fractures.begin();
	while (currentFracture != fractures.end()) {	
		currentFracture->calculate();
		currentFracture++;	
	}
	info("Calculation is done.");
}

Field Stratum::calculateImpactInPoint(const double& x, const double& y) const{
	auto fracture = fractures.begin();
	Field field;
	while (fracture != currentFracture) {
		field += fracture->calculateImpactInPoint(x, y);
		fracture++;
	}
	field.Sxx += Sxx;
	field.Sxy += Sxy;
	field.Syy += Syy;
	return field;
}

void Stratum::visualize() const{
	info("Starting visualisation ...");
	mglGraph gr = mglGraph(0, 1200, 800);
	gr.SetRanges(Xmin, Xmax, Ymin, Ymax);
	gr.Axis();

	drawFractures(gr);
	drawField(gr);
//	drawStressDirections(gr);
	
	gr.WriteFrame("fractures.png");
	info("Visualisation is done.");
}

void Stratum::drawDisplacements() const{
	info("Drawing displacements ...");
	auto fracture = fractures.begin();
	int N = fracture->getNumOfLmnts();
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
	breakerOfFirstFracture.getType(_a, b, c, pressureType);
	// analytical solution
	if ( pressureType == "const" ) {
		for (int i = 0; i < N; i++) {
			double l = fracture->getHalfLength();
			double x = Cx.a[i];
			v.a[i] = (Syy - c) / G * (1 - nu) * sqrt(l*l - x*x);
		}
	} else if ( pressureType == "lag" ) {
		for (int i = 0; i < N; i++) {
			double l = fracture->getHalfLength();
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
		int N = fracture->getNumOfLmnts() + 1;
		double *_x = new double[N];
		double *_y = new double[N];
		fracture->getPointsForPlot(_x, _y);
		mglData x;
		mglData y;
		x.Set(_x, N);
		y.Set(_y, N);
		gr.Plot(x, y, ".k");
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
