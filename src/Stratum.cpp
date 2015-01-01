#include "Stratum.hpp"

Stratum::Stratum() {
	info("Stratum is created.");
}

Stratum::~Stratum() {
}

void Stratum::addFracture(int number, double x, double y, double beta, double h_length,
		int numOfBreaks, double a, double b, double c, std::string pressureType,
		std::string tip, std::string rotation) {
	info("Adding fracture number", number, "at (", x, ",", y, ") ...");
	fractures.push_back(Fracture(this, number, h_length,
					numOfBreaks, a, b, c, pressureType, tip, rotation));
	fractures.back().allocateBreaks(x, y, beta);
}

void Stratum::setRheology(double _G, double _nu) {
	G = _G;		nu = _nu;
}

void Stratum::getRheology(double &_G, double &_nu) {
	_G = G;		_nu = nu;
}

void Stratum::setStresses(double _Sxx, double _Sxy, double _Syy) {
	Sxx = _Sxx;
	Sxy = _Sxy;
	Syy = _Syy;
}

void Stratum::setRanges(double _Xmin, double _Xmax, double _Ymin, double _Ymax) {
	Xmin = _Xmin;
	Xmax = _Xmax;
	Ymin = _Ymin;
	Ymax = _Ymax;		
}

void Stratum::reserve(int numberOfFractures) {
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

Field Stratum::calculateImpactInPoint(const double& x, const double& y) {
	std::vector<Fracture>::iterator fracture = fractures.begin();
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

void Stratum::visualize() {
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

void Stratum::drawFractures(mglGraph& gr) {
	info("Drawing fractures ...");
	std::vector<Fracture>::const_iterator fracture = fractures.begin();
	while (fracture != fractures.end()) {
		int N = fracture->getNumOfBreaks() + 1;
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

void Stratum::drawField(mglGraph &gr) {
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

void Stratum::drawStressDirections(mglGraph &gr) {
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

void Stratum::drawDisplacements() {
	info("Drawing displacements ...");
	std::vector<Fracture>::const_iterator fracture = fractures.begin();
	int N = fracture->getNumOfBreaks();
	double *_x = new double[N];
	double *_v = new double[N];
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
	for (int i = 0; i < N; i++) {
		v.a[i] = Syy / G * (1 - nu) * 
			sqrt(fracture->getLength()*fracture->getLength() - Cx.a[i]*Cx.a[i]);
	}
	gr.Plot(Cx, v, "r");
	gr.WriteFrame("displacements.png");
	delete [] _x;
	delete [] _v;
}