#ifndef SPH_H
#define SPH_H
#include "particle.h"
#include "print_particules.h"
#include "kernel.h"
#include "derivatives.h"

typedef struct Setup Setup;
typedef struct Residual Residual;
typedef void(*update_positions)(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup);


struct Setup{
	int itermax;
	double timestep;
	double kh;
	Verlet* verlet;
	Kernel kernel;
};

struct Residual {
	double mass_eq;
	double momentum_x_eq;
	double momentum_y_eq;
};

Setup* Setup_new(int iter, double timestep,double kh, Verlet* verlet, Kernel kernel);
void Setup_free(Setup* setup);

Residual* Residual_new();
void free_Residuals(Residual** residuals,int N);

void simulate(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, update_positions update_positions, Setup* setup, Animation* animation);
void random_moves(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup);
void update_positions_seminar_5(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup);
void update_positions_ellipse(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup);
void update_positions_test_static_bubble(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup);

void compute_Cs(Particle *particle, Kernel kernel, double kh);
void assemble_residual_NS(Particle* particle, Particle_derivatives* Particle_derivatives, Residual* residual);
void time_integrate(Particle* particle, Residual* residual, double delta_t);

void compute_XSPH_correction(Particle *particle, Kernel kernel, double kh);

void assemble_residual_NS_test(Particle* particle, Particle_derivatives* particle_derivatives, Residual* residual, double radius_circle);

double compute_admissible_dt(double safety_param, double h_p, double c_0, double rho_0, double mu, double sigma);

void compute_normal(Particle *particle, Particle_derivatives* particle_derivatives);

#endif
