#include "kernel.h"

 xy* grad_kernel(xy* p1, xy* p2,double kh, Kernel kernel)
{
	double d = sqrt(pow(p1->x - p2->x, 2) + pow(p1->y - p2->y, 2));
	// double d_x = fabs(p1->x-p2->x); // wrong! grad_kernel is not symmetric
	// double d_y = fabs(p1->y-p2->y);
    double d_x = p1->x-p2->x;
    double d_y = p1->y-p2->y;

	double g;
	double h;
	if (kernel == Cubic) {
		h = kh / 2;
		g = derivative_Cubic_kernel(d / h, h);
	}
	else if (kernel == Lucy) {
		h = kh;
		g = derivative_Lucy_kernel(d / h, h);
	}
	else if (kernel == NewQuartic) {
		h = kh / 2;
		g = derivative_NewQuartic_kernel(d / h, h);
	}
	else if (kernel == Quintic) {
		h = kh / 3;
		g = derivative_Quintic_kernel(d / h, h);
	}

	double g_x = g*(d_x / (h*d));
	double g_y = g*(d_y / (h*d));

	return xy_new(g_x, g_y);
}
 double derivative_Cubic_kernel(double q,double h)
 {
	 double alpha = 15.0/(7.0*M_PI*pow(h, 2));
	 double g;
	 if (q >= 0 && q <= 1)
		 g = -2 * q + 1.5*pow(q, 2);
	 else if (q > 1 && q <= 2)
		 g = -0.5*pow((2 - q), 2);
	 else
		 g = 0;
	 return alpha*g;
 }

 double derivative_Lucy_kernel(double q, double h)
 {
	 double alpha = 5.0/(M_PI*pow(h, 2));
	 double g;
	 if (q >= 0 && q <= 1)
		 g = 3 * (1 - pow(q, 2) - 4 * pow(q, 3));
	 else
		 g = 0;
	 return alpha*g;
 }

 double derivative_NewQuartic_kernel(double q, double h)
 {
	 double alpha = 15.0/(7.0*M_PI*pow(h,2));
	 double g;
	 if (q >= 0 && q <= 2)
		 g = (-9.0/4.0)*q+(19.0/8.0)*pow(q,2)-(5.0/8.0)*pow(q,3);
	 else
		 g = 0;
	 return alpha*g;
 }

 double derivative_Quintic_kernel(double q, double h)
 {
	 double alpha = 7.0 / (478.0*M_PI*pow(h, 2));
	 double g;
	 if (q >= 0 && q <= 1)
		 g = -5*pow(3-q,4)+30*pow(2-q,4)-75*pow(1-q,4);
	 else if (q > 1 && q <= 2)
		 g = -5 * pow(3 - q, 4) + 30 * pow(2 - q, 4);
	 else if (q > 2 && q <= 3)
		 g = -5 * pow(3 - q, 4);
	 else
		 g = 0;
	 return alpha*g;
 }
