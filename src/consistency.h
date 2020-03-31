#ifndef CONSISTENCY_H
#define CONSISTENCY_H
#include "particle.h"
#include "print_particules.h"
#include "kernel.h"
#include "derivatives.h"
#include "SPH.h"

xy* correct_grad(xy *current_grad, Particle *p, double kh, Kernel kernel);
void density_correction_MLS(Particle* pi, double kh, Kernel kernel);
double get_W_MLS(Particle* pi, Particle* pj, double kh, Kernel kernel, double* beta);
double** get_A(Particle* pi, double kh, Kernel kernel);
double* get_beta(double** A);
void Corrective_Smoothed_Particle_Method(Particle *p,Particle_derivatives *dp, double kh, Kernel kernel);
#endif
