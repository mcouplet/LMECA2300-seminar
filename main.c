#include "print_particules.h"
#include "Particule.h"
#include "SPH.h"
#include <math.h>
#include "crtdbg.h" //for memory leak detection

//sans verlet
void script1()
{
	int N = 2000;
	Particle* particles = build_Particles(N, 100);
	double kh = 30;
	double timestep = 1;
	Verlet* verlet = NULL;
	Grid* grid = Grid_new(-100, 100, -100, 100, kh);
	Animation* animation = Animation_new(N, 0.2, grid);
	Kernel kernel = Cubic;
	Setup* setup = Setup_new(50, timestep, verlet,kernel);

	simulate(grid, particles, N, setup, animation);

	Particle_free(particles, N);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);
}

void script2()
{
	int N = 2000;
	Particle* particles = build_Particles(N, 100);
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

	Particle_free(particles, N);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);
}
int main()
{
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	script1();

	return EXIT_SUCCESS;
}