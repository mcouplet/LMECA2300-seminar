#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include "utils.h"
#include "particle.h"
#include "kernel.h"

typedef void*(*getter)(const Particle* particle); // getter function
typedef xy*(*grad_W)(const Particle* particle); // getter function

double compute_div(Particle * particle, getter get, Kernel kernel, double kh);

#endif
