#ifndef CONSISTENCY_H
#define CONSISTENCY_H
#include "particle.h"
#include "print_particules.h"
#include "kernel.h"
#include "derivatives.h"
#include "SPH.h"

//double density_correction_MLS(Particle** particles, Setup* setup, int n_p);
xy* correct_grad(xy *current_grad, Particle *p, double kh, Kernel kernel);

#endif
