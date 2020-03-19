#include "SPH.h"

Setup* Setup_new(int iter, double timestep,Verlet* verlet,Kernel kernel)
{
	Setup* setup = (Setup*)malloc(sizeof(Setup));
	setup->itermax = iter;
	setup->timestep = timestep;
	setup->verlet = verlet;
	setup->kernel = kernel;
	return setup;
}

void Setup_free(Setup* setup)
{
	if (setup->verlet != NULL)
		free(setup->verlet);
	free(setup);
}

void simulate(Grid* grid, Particle* particles, int N, Setup* setup,Animation* animation)
{
	double maxspeed = 2;
	for (int iter = 0;iter < setup->itermax;iter++)
	{
		updateCells(grid, particles, N);
		update_neighborhoods(grid, particles, N, iter,setup->verlet);

		if (animation != NULL)
			display_particles(particles, animation, false);
		random_moves(grid, particles, N, setup->timestep, maxspeed);

	}
	updateCells(grid, particles, N);
	update_neighborhoods(grid, particles, N, 0, setup->verlet);
	if (animation != NULL)
		display_particles(particles, animation, true);
}



//move randomly each particles
void random_moves(Grid* grid, Particle* particles, int N, double timestep, double maxspeed)
{
	for (int i = 0; i < N; i++) {
		double angle = rand_interval(0, 2)*M_PI;
		double speed = rand_interval(0, maxspeed);
		Particle p = particles[i];
		p.v->x = speed * cos(angle);
		p.v->y = speed * sin(angle);
		p.pos->x += p.v->x * timestep;
		p.pos->y += p.v->y * timestep;

		double s = 2;
		//bouncing with the wall
		if (p.pos->x < grid->left)
			p.pos->x = grid->left + s;
		if (p.pos->x > grid->right)
			p.pos->x = grid->right - s;
		if (p.pos->y < grid->bottom)
			p.pos->y = grid->bottom + s;
		if (p.pos->y > grid->top)
			p.pos->y = grid->top - s;
	}
}