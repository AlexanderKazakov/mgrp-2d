#include "Visualization.hpp"

Visualization::Visualization() {
	gr.SetRanges(-2, 20, -20, 20);
	gr.Axis();
}

Visualization::~Visualization() {
	gr.WriteFrame("fractures.png");
}

void Visualization::plotFracture(int N, double* _x, double* _y) {
	mglData x;
	mglData y;
	x.Set(_x, N);
	y.Set(_y, N);
	
	gr.Plot(x, y, "k");
}


void Visualization::sample1() {
	mglData dat(30, 40);	// data to for plotting
	for (long i = 0; i < 30; i++)
		for (long j = 0; j < 40; j++)
			dat.a[i + 30 * j] = 1 / (1 + (i - 15)*(i - 15) / 225. + (j - 20)*(j - 20) / 400.);
	mglGraph gr;	// class for plot drawing
	gr.Rotate(50, 60);	// rotate axis
	gr.Light(true);	// enable lighting
	gr.Surf(dat);	// plot surface
	gr.Cont(dat,"y");	// plot yellow contour lines
	gr.Axis();	// draw axis
	gr.WriteFrame("sample1.png");	// save it
}

void Visualization::sample2() {
	mglData dat(20, 20);	// data to for plotting
	for (long i = 0; i < 20; i++)
		for (long j = 0; j < 20; j++)
			dat.a[i + 20 * j] = 20*sin(i*j);
	mglData x(20);
	mglData y(20);
	mglData z(20);
	mglData s(20);
	for (long i = 0; i < 20; i++) {
			x.a[i] = i;
			s.a[i] = 10;
			y.a[i] = i*i;
			z.a[i] = i / 3.;
	}
	mglGraph gr;	// class for plot drawing
	gr.SetRanges(0, 20, 0, 20);
	gr.Dens(x, x, dat);	// plot surface
	//gr.Cont(dat,"y");	// plot yellow contour lines
	gr.Plot(x, y, "k");
	gr.Plot(x, z, "k");
	gr.Plot(s, x, "k");
	gr.Axis();	// draw axis
	gr.Colorbar();
	gr.WriteFrame("sample2.png");	// save it
}

