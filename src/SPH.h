#ifndef SPH_H
#define SPH_H
#include "particle.h"
#include "print_particules.h"
#include "kernel.h"

typedef struct Setup Setup;

struct Setup{
	int itermax;
	double timestep;
	Verlet* verlet;
	Kernel kernel;
};

Setup* Setup_new(int iter, double timestep, Verlet* verlet, Kernel kernel);
void Setup_free(Setup* setup);

void simulate(Grid* grid, Particle** particles, int N, Setup* setup, Animation* animation);

void random_moves(Grid* grid, Particle** particles, int N, double timestep, double maxspeed);

void assemble_residual_NS(Particle* particle, Particle_derivatives* Particle_derivatives, Residual* residual);

void time_integrate(Particle* particle, Residual* residual, double delta_t);

#endif
