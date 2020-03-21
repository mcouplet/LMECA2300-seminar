#ifndef KERNEL_H
#define KERNEL_H
#include "utils.h"

typedef enum Kernel Kernel;

enum Kernel {Cubic,Lucy,NewQuartic,Quintic};

xy* grad_kernel(xy* p1, xy* p2, double kh, Kernel kernel);
double eval_kernel(xy *p1, xy *p2, double kh, Kernel kernel);

// TODO: make these private
double derivative_Cubic_kernel(double q, double h);
double derivative_Lucy_kernel(double q, double h);
double derivative_NewQuartic_kernel(double q, double h);
double derivative_Quintic_kernel(double q, double h);
#endif
