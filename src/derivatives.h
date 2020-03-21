#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include "utils.h"
#include "particle.h"
#include "kernel.h"

typedef double(*scalar_getter)(Particle* particle); // scalar getter function (e.g. pressure getter)
typedef xy*(*xy_getter)(Particle* particle); // xy getter function (e.g. velocity getter)

//typedef void*(*getter)(const Particle* particle); // getter function
//typedef xy*(*grad_W)(const Particle* particle); // getter function

double compute_div(Particle * particle, xy_getter get, Kernel kernel, double kh);
xy * compute_grad(Particle * particle, scalar_getter get, Kernel kernel, double kh);
double compute_lapl(Particle *particle, scalar_getter get, Kernel kernel, double kh);

#endif
