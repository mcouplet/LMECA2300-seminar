#include "print_particules.h"
#include "particle.h"
#include "SPH.h"
#include "derivatives.h"
#include <math.h>
#include "kernel.h"
#include "consistency.h"

//#include "crtdbg.h" // for memory leak detection; comment if you're on Linux

void script_csf();
void script_csf_circle();
void script_circle_to_ellipse();
void script_csf_circle_paper();
void script1();
void script2();

int main() {
	// _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF); // comment if on Linux

	dam_break();

	return EXIT_SUCCESS;
}
void dam_break() {

	// Parameters of the problem
	double L = 1; // size of the domain: [-L,L] x [-L,L]
	double l = 0.4; // size of the square: [-l,l] x [-l,l]
	double dt = 0.001; // physical time step
	double T = 30; // duration of simulation
	// double T = dt; // duration of simulation

	// Physical parameters
	double rho_0 = 998.0; // initial (physical) density of water at 20°C (in kg/m^3)
	double mu = 1.0016e-3; // dynamic viscosity of water at 20°C (in N.s/m^2)
	double gamma = 7; // typical value for liquid (dimensionless)
	double c_0 = 1.0;//1481; // sound speed in water at 20°C (in m/s)
	double sigma = 72.86e-3; // surface tension of water-air interface at 20°C (in N/m)
	bool gravity = 1;

	// SPH parameters
	int n_per_dim = 51; // number of particles per dimension
	double kh = sqrt(21) * 2 * l / n_per_dim; // kernel width to ensure 21 particles in the neighborhood
	int n_iter = (int)(T / dt); // number of iterations to perform
	Kernel kernel = Cubic; // kernel choice
	double interface_threshold = 1.5;//20; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)
	Verlet *verlet = NULL; // don't use Verlet (for now)
	double XSPH_epsilon = 0.5;
	Free_surface_detection surface_detection = DIVERGENCE;
	double CR = 0.5;
	double CF = 0.5;
	// Free_surface_detection surface_detection = CSF;


	printf("n_iter = %d\n", n_iter);

	// Animation parameter
	double T_anim = 10; // duration of animation
	double dt_anim = T_anim / n_iter; // time step of animation

	// Initialize particles on a square
	int n_p = squared(n_per_dim); // total number of particles
	double h = 2 * l / (n_per_dim - 1); // step between neighboring particles
	double m = rho_0 * h*h;

	Particle** particles = (Particle**)malloc(n_p * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(n_p * sizeof(Particle_derivatives*));
	Residual** residuals = malloc(n_p * sizeof(Residual*));
	for (int i = 0; i < n_per_dim; i++) {
		for (int j = 0; j < n_per_dim; j++) {
			int index = i * n_per_dim + j;
			xy *pos = xy_new(-l + i * h, -l + j * h);
			xy *v = xy_new(0, 0); // initial velocity = 0
			particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
		}
	}

	// Setup grid
	Grid *grid = Grid_new(-L, L, -L, L, kh);
	// Setup BOUNDARY
	double d = 0.5*L;
	Boundary* boundary = Boundary_new(-l,d,-l,d,CR,CF);
	// Setup setup
	Setup *setup = Setup_new_bis(n_iter, dt, kh, verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon, gravity);
	// Setup animation
	Animation *animation = Animation_new(n_p, dt_anim, grid, 1);

	// Start simulation
	simulate_boundary(grid, particles, particles_derivatives, residuals, n_p, update_positions_seminar_5, setup, animation, boundary);
	// simulate(grid, particles, particles_derivatives, residuals, n_p, update_positions_seminar_5, setup, animation);
	// Free stuff
	free_particles(particles, n_p);
	free_particles_derivatives(particles_derivatives, n_p);
	free_Residuals(residuals, n_p);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);

}
