#ifndef SPH_H
#define SPH_H
#include "particle.h"
#include "print_particules.h"
#include "kernel.h"

typedef struct Setup Setup;
typedef enum Free_surface_detection Free_surface_detection;

struct Setup{
	int itermax;
	double timestep;
	Verlet* verlet;
	Kernel kernel;
};

enum Free_surface_detection {
    CSF = 1,
    DIVERGENCE = 2
};

Setup* Setup_new(int iter, double timestep, Verlet* verlet, Kernel kernel);
void Setup_free(Setup* setup);

void simulate(Grid* grid, Particle** particles, int N, Setup* setup, Animation* animation);

void random_moves(Grid* grid, Particle** particles, int N, double timestep, double maxspeed);

void assemble_residual_NS(Particle* particle, Particle_derivatives* particle_derivatives, Residual* residual, Free_surface_detection detection);

void time_integrate(Particle* particle, Residual* residual, double delta_t);

void compute_Cs(Particle *particle, Kernel kernel, double kh);

void compute_XSPH_correction(Particle *particle, Kernel kernel, double kh);

void assemble_residual_NS_test(Particle* particle, Particle_derivatives* particle_derivatives, Residual* residual, double radius_circle);

#endif
