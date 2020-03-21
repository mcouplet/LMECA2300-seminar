#include "SPH.h"

Setup* Setup_new(int iter, double timestep,Verlet* verlet,Kernel kernel) {
	Setup* setup = (Setup*)malloc(sizeof(Setup));
	setup->itermax = iter;
	setup->timestep = timestep;
	setup->verlet = verlet;
	setup->kernel = kernel;
	return setup;
}

void Setup_free(Setup* setup) {
	if(setup->verlet != NULL)
		free(setup->verlet);
	free(setup);
}

void simulate(Grid* grid, Particle** particles, int N, Setup* setup, Animation* animation) {
	double maxspeed = 2;
	for(int iter = 0; iter < setup->itermax; iter++) {
		update_cells(grid, particles, N);
		update_neighborhoods(grid, particles, N, iter, setup->verlet);

		if (animation != NULL)
			display_particles(particles, animation, false);

		

		//random_moves(grid, particles, N, setup->timestep, maxspeed); // move particles randomly
	}
	update_cells(grid, particles, N);
	update_neighborhoods(grid, particles, N, 0, setup->verlet);
	if (animation != NULL)
		display_particles(particles, animation, true);
}



//move randomly each particles
void random_moves(Grid* grid, Particle** particles, int N, double timestep, double maxspeed) {
	for (int i = 0; i < N; i++) {
		double angle = rand_interval(0, 2)*M_PI;
		double speed = rand_interval(0, maxspeed);
		Particle *p = particles[i];
		p->v->x = speed * cos(angle);
		p->v->y = speed * sin(angle);
		p->pos->x += p->v->x * timestep;
		p->pos->y += p->v->y * timestep;

		double s = 2;
		//bouncing with the wall
		if (p->pos->x < grid->left)
			p->pos->x = grid->left + s;
		if (p->pos->x > grid->right)
			p->pos->x = grid->right - s;
		if (p->pos->y < grid->bottom)
			p->pos->y = grid->bottom + s;
		if (p->pos->y > grid->top)
			p->pos->y = grid->top - s;
	}
}


// Assemble the residual of the (incompressible) Navier-Stokes equations based on the derivatives available
void assemble_residual_NS(Particle* particle, Particle_derivatives* particle_derivatives, Residual* residual) {
  double mu_i = particle->param->dynamic_viscosity;
  
  double rho_i = particle->rho;
  double div_vel_i = particle_derivatives->div_v;
  xy* grad_P = particle_derivatives->grad_P;
  xy* lapl_v = particle_derivatives->lapl_v;
  
  xy* f_CSF = xy_new(0.0, 0.0); // TODO: compute_CSF_term(particle, particle_derivatives);
  
  residual->mass_eq = -rho_i * div_vel_i;
  residual->momentum_x_eq = (-1.0/rho_i) * grad_P->x + (mu_i/rho_i) * lapl_v->x + f_CSF->x;
  residual->momentum_y_eq = (-1.0/rho_i) * grad_P->y + (mu_i/rho_i) * lapl_v->y + f_CSF->y;
  
}

// Time integrate the Navier-Stokes equations based on the residual already assembled
void time_integrate(Particle* particle, Residual* residual, double delta_t) {
 
  particle->rho += delta_t * residual->mass_eq;
  particle->v->x += delta_t * residual->momentum_x_eq;
  particle->v->y += delta_t * residual->momentum_y_eq;
  
}







