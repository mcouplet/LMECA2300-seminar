#include "print_particules.h"
#include "particle.h"
#include "SPH.h"
#include "derivatives.h"
#include <math.h>
#include "kernel.h"
//#include "crtdbg.h" // for memory leak detection; comment if you're on Linux

void script_csf();
void script_csf_circle();
void script_circle_to_ellipse();
void script_csf_circle_paper();
void script1();
void script2();

int main() {
	// _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF); // comment if on Linux
	script_csf();
	//script_csf_circle();
	//script_circle_to_ellipse();
// 	script_csf_circle_paper();
	//script1();

	return EXIT_SUCCESS;
}

// Evolution of a 2D square submitted to the tension surface force
void script_csf() {

	// Parameters of the problem
	double l1 = 0.01*sqrt(M_PI)/2*10; // extent of ethanol: [-l,l] x [-l,l]
	double l0 = 1.5*l1; // extent of air: [-L,L] x [-L,L]
	double dt = 0.001; // physical time step
	double T = 1; // duration of simulation
	// double T = dt; // duration of simulation

	// Physical parameters
	double rho0 = 1.117; // air density
	double rho1 = 785.92; // ethanol density
	double mu0 = 1.81e-5; // air viscosity
	double mu1 = 1.261e-3; // ethanol viscosity),
	double sigma = 22.39e-3; // (surface tension),
	double gamma = 1; // typical value for liquid (dimensionless)
	double c0 = 343; // sound speed in air
	double c1 = 1144; // sound speed in ethanol

	double L = 1.1*l0; // size of grid (allow some margin)

	// SPH parameters
	int n_per_dim = 151; // number of particles per dimension
	double kh = sqrt(21) * 2*L / n_per_dim; // kernel width to ensure 21 particles in the neighborhood
	int n_iter = (int)(T / dt); // number of iterations to perform
	Kernel kernel = Cubic; // kernel choice
	Verlet *verlet = NULL; // don't use Verlet (for now)
	double XSPH_epsilon = 0.5;
	Free_surface_detection surface_detection = CSF;
	double interface_threshold = 20; //20; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)

	printf("n_iter = %d\n", n_iter);
	printf("l0 = %lf, l1 = %lf\n", l0, l1);


	// Animation parameter
	// double T_anim = 10000*T; // duration of animation
	// double dt_anim = T_anim / n_iter; // time step of animation
	double dt_anim = 0.1;

	// Initialize particles on a square
	int n_p = squared(n_per_dim); // total number of particles
	double h = 2 * l0 / (n_per_dim - 1); // step between neighboring particles
	double m0 = rho0 * h*h; // mass per particle of air
	double m1 = rho1 * h*h; // mass per particle of ethanol

	double k_max = M_PI / h;
	double c_phi = sqrt(sigma * k_max / (rho0 + rho1));
	if(c_phi * dt / h < 0.5) printf("Stability condition is verified.\n");
	else printf("WARNING: Stability condition is not verified!\n");

	Particle** particles = (Particle**)malloc(n_p * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(n_p * sizeof(Particle_derivatives*));
	Residual** residuals = malloc(n_p * sizeof(Residual*));
	for (int i = 0; i < n_per_dim; i++) {
		for (int j = 0; j < n_per_dim; j++) {
			int index = i * n_per_dim + j;
			xy *pos = xy_new(-l0 + i * h, -l0 + j * h);
			xy *v = xy_new(0, 0); // initial velocity = 0
			if(fabs(pos->x) <= l1 && fabs(pos->y) <= l1) // ethanol
				particles[index] = Particle_new(index, m1, pos, v, rho1, mu1, c1, gamma, sigma);
			else // air
				particles[index] = Particle_new(index, m0, pos, v, rho0, mu0, c0, gamma, sigma);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
		}
	}

	// Setup grid
	Grid *grid = Grid_new(-L, L, -L, L, kh);
	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon);
	// Setup animation
	Animation *animation = Animation_new(n_p, dt_anim, grid, 1);

	// Start simulation
	simulate(grid, particles, particles_derivatives, residuals, n_p, update_positions_seminar_5, setup, animation);

	// Free stuff
	free_particles(particles, n_p);
	free_particles_derivatives(particles_derivatives, n_p);
	free_Residuals(residuals, n_p);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);
}

// // Equilibrium rod
// void script_rod() {
// 	// Parameters of the problem
// 	double r = 0.25; // radius of disk
// 	double l = 1; // width of square domain
// 	double dt = 0.001; // physical time step
// 	// double T = 1; // duration of simulation
// 	double T = dt; // duration of simulation
//
// 	// Physical parameters
// 	double rho0 = 1.117; // air density
// 	double rho1 = 785.92; // ethanol density
// 	double mu0 = 1.81e-5; // air viscosity
// 	double mu1 = 1.261e-3; // ethanol viscosity),
// 	double sigma = 22.39e-3; // (surface tension),
// 	double gamma = 1; // typical value for liquid (dimensionless)
// 	double c0 = 343; // sound speed in air
// 	double c1 = 1144; // sound speed in ethanol
//
// 	double L = 1.1*l0; // size of grid (allow some margin)
//
// 	// SPH parameters
// 	int n_per_dim = 151; // number of particles per dimension
// 	double kh = sqrt(21) * 2*L / n_per_dim; // kernel width to ensure 21 particles in the neighborhood
// 	int n_iter = (int)(T / dt); // number of iterations to perform
// 	Kernel kernel = Cubic; // kernel choice
// 	Verlet *verlet = NULL; // don't use Verlet (for now)
// 	double XSPH_epsilon = 0.5;
// 	Free_surface_detection surface_detection = CSF;
// 	double interface_threshold = 20; //20; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)
//
// 	printf("n_iter = %d\n", n_iter);
// 	printf("l0 = %lf, l1 = %lf\n", l0, l1);
//
//
// 	// Animation parameter
// 	// double T_anim = 10000*T; // duration of animation
// 	// double dt_anim = T_anim / n_iter; // time step of animation
// 	double dt_anim = 0.1;
//
// 	// Initialize particles on a square
// 	int n_p = squared(n_per_dim); // total number of particles
// 	double h = 2 * l0 / (n_per_dim - 1); // step between neighboring particles
// 	double m0 = rho0 * h*h; // mass per particle of air
// 	double m1 = rho1 * h*h; // mass per particle of ethanol
//
// 	double k_max = M_PI / h;
// 	double c_phi = sqrt(sigma * k_max / (rho0 + rho1));
// 	if(c_phi * dt / h < 0.5) printf("Stability condition is verified.\n");
// 	else printf("WARNING: Stability condition is not verified!\n");
//
// 	Particle** particles = (Particle**)malloc(n_p * sizeof(Particle*));
// 	Particle_derivatives** particles_derivatives = malloc(n_p * sizeof(Particle_derivatives*));
// 	Residual** residuals = malloc(n_p * sizeof(Residual*));
// 	for (int i = 0; i < n_per_dim; i++) {
// 		for (int j = 0; j < n_per_dim; j++) {
// 			int index = i * n_per_dim + j;
// 			xy *pos = xy_new(-l0 + i * h, -l0 + j * h);
// 			xy *v = xy_new(0, 0); // initial velocity = 0
// 			if(fabs(pos->x) <= l1 && fabs(pos->y) <= l1) // ethanol
// 				particles[index] = Particle_new(index, m1, pos, v, rho1, mu1, c1, gamma, sigma);
// 			else // air
// 				particles[index] = Particle_new(index, m0, pos, v, rho0, mu0, c0, gamma, sigma);
// 			particles_derivatives[index] = Particle_derivatives_new(index);
// 			residuals[index] = Residual_new();
// 		}
// 	}
//
// 	// Setup grid
// 	Grid *grid = Grid_new(-L, L, -L, L, kh);
// 	// Setup setup
// 	Setup *setup = Setup_new(n_iter, dt, kh, verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon);
// 	// Setup animation
// 	Animation *animation = Animation_new(n_p, dt_anim, grid, 1);
//
// 	// Start simulation
// 	simulate(grid, particles, particles_derivatives, residuals, n_p, update_positions_seminar_5, setup, animation);
//
// 	// Free stuff
// 	free_particles(particles, n_p);
// 	free_particles_derivatives(particles_derivatives, n_p);
// 	free_Residuals(residuals, n_p);
// 	Grid_free(grid);
// 	Setup_free(setup);
// 	Animation_free(animation);
// }

// Evolution of a 2D circle submitted to the tension surface force
void script_csf_circle() {

	// Computational Parameters
	double L = 100; // size of the domain: [-L,L] x [-L,L]
	double l = 50; // size of the square: [-l,l] x [-l,l]
	int n_per_dim = 51; // number of particles per dimension
	double kh = 0.1*L; // kernel width
	double dt = 1.0; // physical time step
	double dt_anim = 0.001; // time step for animation
	int n_iter = 2000; // number of iterations to perform
	Verlet *verlet = NULL; // don't use Verlet (for now)
	Kernel kernel = Cubic; // kernel choice
	double interface_threshold = 0.08; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)
	double XSPH_epsilon = 0.0; //0.5;
	Free_surface_detection surface_detection = CSF;

	// Physical parameters
	double rho_0 = 998.0; // initial (physical) density of water at 20°C (in kg/m^3)
	double mu = 1.0016e-3; // dynamic viscosity of water at 20°C (in N.s/m^2)
	double gamma = 7.0; // typical value for liquid (dimensionless)
	double c_0 = 1.0;//1481; // sound speed in water at 20°C (in m/s)
	double sigma = 72.86e-3; // surface tension of water-air interface at 20°C (in N/m)

	// Initialize particles in a circle
	int n_p = squared(n_per_dim); // total number of particles
	double h = 2 * l / (n_per_dim - 1); // step between neighboring particles
	double m = rho_0 * M_PI * l * l / n_p;
	double deltat_max = 0.25 * h / c_0;
	// 	printf(">>>>>>>>>>>>>> deltat_max = %2.6f \n", deltat_max);

	Particle** particles = (Particle**)malloc(n_p * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(n_p * sizeof(Particle_derivatives*));
	Residual** residuals = (Residual**)malloc(n_p * sizeof(Residual*));

	int k = 0;
	int alpha = 2;
	int nb = (int)(alpha*n_per_dim);
	for (int i = 0; i < n_per_dim; i++) {
		for (int j = 0; j < n_per_dim; j++) {
			int index = i * n_per_dim + j;
			xy *pos = generate_circle(k, n_p, nb, l);
			// 			xy *pos_circle = map_to_circle(pos);
			// 			pos->x = pos_circle->x;
			// 			pos->y = pos_circle->y;
			xy *v = xy_new(0, 0); // initial velocity = 0
			particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
			k++;
		}
	}


	// Setup grid
	Grid *grid = Grid_new(-L, L, -L, L, kh);
	// Setup animation
	Animation *animation = Animation_new(n_p, dt_anim, grid, 1);
	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon);

	// Start simulation
	simulate(grid, particles, particles_derivatives, residuals, n_p, update_positions_seminar_5, setup, animation);

	// Free stuff
	free_particles(particles, n_p);
	free_particles_derivatives(particles_derivatives, n_p);
	free_Residuals(residuals, n_p);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);

}

// Evolution of a 2D circle with non-zero initial velocities (no surface tension force)
// Test case from "Simulating Free Surface Flows with SPH", Monaghan (1994)
void script_circle_to_ellipse() {

	// Parameters of the problem
	double l = 1.0; // radius of the circle
	double L = 2.0*l; // size of the domain: [-L,L] x [-L,L]
	double dt = 1.0e-5; // physical time step
	double T = 0.0076; // duration of simulation

	// Physical parameters
	double rho_0 = 1000.0; // initial (physical) density of water at 20°C (in kg/m^3)
	double mu = 1.0016e-3; // dynamic viscosity of water at 20°C (in N.s/m^2)
	double gamma = 7.0; // typical value for liquid (dimensionless)
	double c_0 = 1400.0;//1481; // sound speed in water at 20°C (in m/s)
	double sigma = 0.0; // surface tension of water-air interface at 20°C (in N/m)

	// SPH parameters
	Verlet *verlet = NULL; // don't use Verlet (for now)
	Kernel kernel = Cubic; // kernel choice
	double interface_threshold = 1.5;//1000.0; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)
	double XSPH_epsilon = 0.5;
	Free_surface_detection surface_detection = DIVERGENCE;
	int N_c = 30; // number of circonferences on which points are placed
	int N_p = 6; // number of points on the first circonference (doubled for every circonference)
	int N_tot = 1; // total number of points
	for (int i = 1; i < N_c; i++) {
		N_tot += i * N_p;
	}
	printf("N_tot = %d \n", N_tot);
	int n_iter = (int)(T / dt); // number of iterations to perform
	double kh = 0.2*l;// is ideal to reach t = 0.0076; // kernel width


	// Animation parameter
	double T_anim = 0.1; // duration of animation
	double dt_anim = T_anim / n_iter; // time step of animation

	// Initialize particles in a circle
	double m = rho_0 * M_PI * l * l / N_tot; // mass of each particle

	Particle** particles = (Particle**)malloc(N_tot * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(N_tot * sizeof(Particle_derivatives*));
	Residual** residuals = malloc(N_tot * sizeof(Residual*));

	// parameters defining the circle
	double b, delta_s, k, theta;
	delta_s = l / ((double)N_c - 1.0);
	theta = (2*M_PI) / ((double)N_p);

	int index = 0;
	for (int i = 0; i < N_c; i++) {
		b = i;
		if (b == 0) {
			xy *pos = xy_new(0.0, 0.0);
			xy *v = xy_new(0.0, 0.0);
			particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
			index++;
		}
		else {
			for (int j = 0; j < i*N_p; j++) {
				k = (double)j / b;
				xy *pos = xy_new(b*delta_s*cos(k*theta), b*delta_s*sin(k*theta));
				xy *v = xy_new(-100.0*pos->x, 100.0*pos->y);
				particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma);
				particles_derivatives[index] = Particle_derivatives_new(index);
				residuals[index] = Residual_new();
				index++;
			}
		}
	}


	// Setup grid
	Grid *grid = Grid_new(-L, L, -L, L, kh);
	// Setup animation
	Animation *animation = Animation_new(N_tot, dt_anim, grid, 1);
	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon);

	// Start simulation
	simulate(grid, particles, particles_derivatives, residuals, N_tot, update_positions_ellipse, setup, animation);

	// Free stuff
	free_particles(particles, N_tot);
	free_particles_derivatives(particles_derivatives, N_tot);
	free_Residuals(residuals, N_tot);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);
}

// Evolution of a 2D circle : test case based on "Modelling surface tension of two-dimensional droplet using
// Smoothed Particle Hydrodynamics", Nowoghomwenma, 2018

void script_csf_circle_paper() {

//	// Nowoghomwenma's paper
// 	// Parameters of the problem
// 	double l = 1e-3; // radius of the circle
// 	double L = 2.0*l; // size of the domain: [-L,L] x [-L,L]
// 	double dt = 0.000001; // physical time step
// 	double T = 0.01; // duration of simulation
//
// 	// Physical parameters
// 	double rho_0 = 1000.0; // initial (physical) density of water at 20°C (in kg/m^3)
// 	double mu = 0.01; // dynamic viscosity of water at 20°C (in N.s/m^2)
// 	double gamma = 7.0; // typical value for liquid (dimensionless)
// 	double c_0 = 5.0;//5.0;//1481; // sound speed in water at 20°C (in m/s)
// 	double sigma = 1.0; // surface tension of water-air interface at 20°C (in N/m)

	// Brackbill's paper
  	// Parameters of the problem
	double l = 2e-2; // radius of the circle
	double L = 2.0*l; // size of the domain: [-L,L] x [-L,L]
	double dt = 0.001; // physical time step
	double T = 1.0; // duration of simulation

	// Physical parameters
	double rho_0 = 1000.0; // initial (physical) density of water at 20°C (in kg/m^3)
	double mu = 0.01;//0.01; // dynamic viscosity of water at 20°C (in N.s/m^2)
	double gamma = 7.0; // typical value for liquid (dimensionless)
	double c_0 = 10.0;//5.0;//1481; // sound speed in water at 20°C (in m/s)
	double sigma = 23.61e-3; // surface tension of water-air interface at 20°C (in N/m)

	// SPH parameters
	Verlet *verlet = NULL; // don't use Verlet (for now)
	Kernel kernel = Cubic; // kernel choice
	double interface_threshold = 1.0; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)
	double XSPH_epsilon = 1.0;
	Free_surface_detection surface_detection = DIVERGENCE;
	int N_c = 11;//21; // number of circonferences on which points are placed
	int N_p = 6; // number of points on the first circonference (doubled for every circonference)
	int N_tot = 1; // total number of points
	for (int i = 1; i < N_c; i++) {
		N_tot += i * N_p;
	}
	printf("N_tot = %d \n", N_tot);
	int n_iter = (int)(T / dt); // number of iterations to perform
	double kh = 1.0*l; // kernel width

	// Animation parameter
	double T_anim = 10; // duration of animation
	double dt_anim = T_anim / n_iter; // time step of animation


	// Initialize particles in a circle
	double m = rho_0 * M_PI * l * l / N_tot; // mass of each particle

	Particle** particles = (Particle**)malloc(N_tot * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(N_tot * sizeof(Particle_derivatives*));
	Residual** residuals = malloc(N_tot * sizeof(Residual*));

	// parameters defining the circle
	double b, delta_s, k, theta;
	delta_s = l / ((double)N_c - 1.0);
	theta = (2*M_PI) / ((double)N_p);
	double ecc = 0.6;

	int index = 0;
	for(int i = 0; i < N_c; i++) {
	  b = i;
	  if (b==0) {
	    xy *pos = xy_new(0.0, 0.0);
	    xy *v = xy_new(0.0, 0.0);
	    particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma);
	    particles_derivatives[index] = Particle_derivatives_new(index);
	    residuals[index] = Residual_new();
	    index++;
	  }
	  else {
	    for (int j = 0; j < i*N_p; j++) {
		k = (double) j/b;
		xy *pos = xy_new(b*delta_s*cos(k*theta), b*delta_s*sin(k*theta));
		if (ecc != 0.0) {
		  pos->x *= sqrt(2.0/sin(ecc*M_PI)) * sin(0.5*ecc*M_PI);
		  pos->y *= sqrt(2.0/sin(ecc*M_PI)) * cos(0.5*ecc*M_PI);
		}
		xy *v = xy_new(0.0, 0.0);
		particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma);
		particles_derivatives[index] = Particle_derivatives_new(index);
		residuals[index] = Residual_new();
// 		if (b == N_c - 1) particles[index]->on_free_surface = true; // WARNING: to be removed
		index++;
	    }
	  }
	}

	//Estimate maximum admissible time step for stability
	double h_p = l / sqrt(N_tot);
	double safety_param = 0.8;
	dt = compute_admissible_dt(safety_param, h_p, c_0, rho_0, mu, sigma);
	n_iter = (int)(T/dt);

	// Setup grid
	Grid *grid = Grid_new(-L, L, -L, L, kh);
	// Setup animation
	Animation *animation = Animation_new(N_tot, dt_anim, grid, 1);
	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon);

	simulate(grid, particles, particles_derivatives, residuals, N_tot, update_positions_seminar_5, setup, animation);
// 	simulate(grid, particles, particles_derivatives, residuals, N_tot, update_positions_test_static_bubble, setup, animation);

	// Free stuff
	free_particles(particles, N_tot);
	free_particles_derivatives(particles_derivatives, N_tot);
	free_Residuals(residuals, N_tot);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);

}

// Without Verlet
void script1() {
	int n_p = 2000;
	Particle** particles = build_particles(n_p, 100);
	double kh = 30;
	double timestep = 1;
	Verlet* verlet = NULL;
	Grid* grid = Grid_new(-100, 100, -100, 100, kh);
	Animation* animation = Animation_new(n_p, 0.05, grid, 1);
	Kernel kernel = Cubic;
	Setup *setup = Setup_new(50, timestep, kh, verlet, kernel, NONE, 0, 0);


	simulate(grid, particles, NULL, NULL, n_p, random_moves, setup, animation);
	free_particles(particles, n_p);
	Setup_free(setup);
	Grid_free(grid);
	Animation_free(animation);
}

// With Verlet
void script2() {
	int n_p = 2000;
	Particle** particles = build_particles(n_p, 100);
	double kh = 10;
	double vmax = 2;
	int T = 4;
	double timestep = 1;
	double L = 2 * T*vmax*timestep;
	Verlet* verlet = Verlet_new(kh, L, T);


	Grid* grid = Grid_new(-100, 100, -100, 100, kh + L);
	Animation* animation = Animation_new(n_p, 0.2, grid, 1);
	Kernel kernel = Cubic;
	Setup *setup = Setup_new(50, timestep, kh, verlet, kernel, NONE, 0, 0);

	simulate(grid, particles, NULL, NULL, n_p, random_moves, setup, animation);
	free_particles(particles, n_p);
	Setup_free(setup);
	Grid_free(grid);
	Animation_free(animation);
}
