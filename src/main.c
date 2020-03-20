#include "print_particules.h"
#include "particle.h"
#include "SPH.h"
#include "derivatives.h"
#include <math.h>
//#include "crtdbg.h" // for memory leak detection; comment if you're on Linux

void script_csf();
void script1();
void script2();

int main() {
	// _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF); // comment if on Linux
	script_csf();

	return EXIT_SUCCESS;
}

// Evolution of a 2D square submitted to the tension surface force
void script_csf() {

	// Parameters
	double L = 100; // size of the domain: [-L,L] x [-L,L]
	double l = 50; // size of the square: [-l,l] x [-l,l]
	int n_per_dim = 51; // number of particles per dimension
	double density = 1; // initial (physical) density
	double Cs = 1; // color field (1 = fluid, 0 = void)
	double kh = 30; // kernel width
	double dt = 1; // physical time step
	double dt_anim = 0.05; // time step for animation
	int n_iter = 50; // number of iterations to perform
	Verlet *verlet = NULL; // don't use Verlet (for now)
	Kernel kernel = Cubic; // kernel choice

	// Initialize particles on a square
	int n_p = squared(n_per_dim); // total number of particles
	double h = 2*l / (n_per_dim-1); // step between neighboring particles
	double m = density * h*h;
	Particle** particles = (Particle**) malloc(n_p*sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(n_p * sizeof(Particle_derivatives*));
	for(int i = 0; i < n_per_dim; i++) {
		for(int j = 0; j < n_per_dim; j++) {
			int index = i*n_per_dim + j;
			xy *pos = xy_new(-l+i*h, -l+j*h);
			xy *v = xy_new(pos->x, pos->y); // initial velocity = 0
			particles[index] = Particle_new(index, m, pos, v, density, Cs);
			particles_derivatives[index] = Particle_derivatives_new(index);
		}
	}

	// Setup grid
	Grid *grid = Grid_new(-L, L, -L, L, kh);
	// Setup animation
	//Animation *animation = Animation_new(n_p, dt_anim, grid);
	Animation *animation = NULL;
	// Setup setup (lol)
	Setup *setup = Setup_new(n_iter, dt, verlet, kernel);


	// Start simulation
	for(int iter = 0; iter < setup->itermax; iter++) {
		update_cells(grid, particles, n_p);
		update_neighborhoods(grid, particles, n_p, iter, setup->verlet);
		if (animation != NULL)
			display_particles(particles, animation, false);

		// Reset derivatives
		for(int i = 0; i < n_p; i++) {
			Particle_derivatives_reset(particles_derivatives[i]);
			particles_derivatives[i]->div_v = compute_div(particles[i], Particle_get_v, kernel, grid->h);
			printf("pos = (%lf,%lf), div_v = %lf\n", particles[i]->pos->x, particles[i]->pos->y, particles_derivatives[i]->div_v);
		}

		// integrate :-)
	}
	update_cells(grid, particles, n_p);
	update_neighborhoods(grid, particles, n_p, 0, setup->verlet);
	if (animation != NULL)
		display_particles(particles, animation, true);

	// Free stuff
	free_particles(particles, n_p);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);

}

// Without Verlet
void script1() {
	int N = 2000;
	Particle** particles = build_particles(N, 100);
	double kh = 30;
	double timestep = 1;
	Verlet* verlet = NULL;
	Grid* grid = Grid_new(-100, 100, -100, 100, kh);
	Animation* animation = Animation_new(N, 0.05, grid);
	Kernel kernel = Cubic;
	Setup* setup = Setup_new(50, timestep, verlet,kernel);

	simulate(grid, particles, N, setup, animation);

}

// With Verlet
void script2() {
	int N = 2000;
	Particle** particles = build_particles(N, 100);
	double kh = 10;
	double vmax = 2;
	int T = 4;
	double timestep = 1;
	double L = 2*T*vmax*timestep;
	Verlet* verlet = Verlet_new(kh, L, T);


	Grid* grid = Grid_new(-100, 100, -100, 100, kh+L);
	Animation* animation = Animation_new(N, 0.2, grid);
	Kernel kernel = Cubic;
	Setup* setup = Setup_new(50, timestep, verlet,kernel);

	simulate(grid, particles, N, setup, animation);

	free_particles(particles, N);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);
}
