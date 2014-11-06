#include "Stratum.hpp"

Stratum::Stratum() {
}

Stratum::~Stratum() {
}

void Stratum::addFracture(const Fracture& fracture) {
	fractures.push_back(fracture);
}

void Stratum::setRheology(double _G, double _nu) {
	G = _G;		nu = _nu;
}

void Stratum::sortFractures() {
	std::sort(fractures.begin(), fractures.end());
	currentFracture = fractures.begin();
}

int Stratum::calculateNextFracture() {
	if (fractures.end() == currentFracture)
		return 0;
	
	(*currentFracture).calculate(fractures.begin());
	currentFracture++;	
	return 1;
}

Field Stratum::calculateImpactInPoint(const double& x, const double& y) {
	std::vector<Fracture>::iterator fracture = fractures.begin();
	Field field;
	while (fracture != currentFracture) {
		field += (*fracture).calculateImpactInPoint(x, y);
		fracture++;
	}
	return field;
}

//void Stratum::visualize() {
//	Visualization vis;
//	std::vector<Fracture>::iterator fracture = fractures.begin();
//	while (fracture != fractures.end()) {
//		int N = (*fracture).getNumOfPointsForPlot();
//		std::cout << N << std::endl;
//		double *x = new double[N];
//		double *y = new double[N];
//
//		(*fracture).getPointsForPlot(x, y);
//		vis.plotFracture(N, x, y);
//		
//		//vis.plotField(&Stratum::calculateImpactInPoint);
//		
//		fracture++;
//		delete [] x;
//		delete [] y;
//	}
//}

void Stratum::visualize() {
	std::vector<Fracture>::iterator fracture = fractures.begin();
	mglGraph gr;
	double Xmin = -2;
	double Xmax = 20;
	double Ymin = -22;
	double Ymax = 22;
	gr.SetRanges(Xmin, Xmax, Ymin, Ymax);
	gr.Axis();
	while (fracture != fractures.end()) {
		int N = (*fracture).getNumOfPointsForPlot();
		std::cout << N << std::endl;
		double *_x = new double[N];
		double *_y = new double[N];

		(*fracture).getPointsForPlot(_x, _y);
		mglData x;
		mglData y;
		x.Set(_x, N);
		y.Set(_y, N);

		gr.Plot(x, y, "k");		
		
		fracture++;
		delete [] _x;
		delete [] _y;
	}
	
//	int N = 100;
//	mglData x(N);
//	mglData y(N);
//	mglData f(N, N);
//	for (int i = 0; i < N; i += 1)
//		for (int j = 0; j < N; j += 1) {
//			x.a[i] = (Xmax - Xmin)*i/N + Xmin;
//			y.a[j] = (Ymax - Ymin)*j/N + Ymin;
//			double _x = x.a[i];
//			double _y = y.a[j];
//			f.a[i + N*j] = 1e+4 * calculateImpactInPoint(_x, _y).Sxy;
//		}
//	gr.Dens(x, y, f);
//	gr.Colorbar();
	gr.WriteFrame("fractures.png");
}