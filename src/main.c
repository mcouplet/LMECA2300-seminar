#include "neighborhood_search.h"
#include "kernel.h"
#include "checkDerivatives.h"
#include <math.h>

int NPTS = 100;

// for v, a floating point value between 0 and 1, this function fills color with
// the improved jet colormap color corresponding to v
static void colormap(float v, float color[3])
{
	float v1 = 3.5 * (v - 0.7);
	float v2 = 1.25 * v;
	float v3 = fminf(0.5, v) * 2.0;

	color[0] = -v1 * v1 + 1.0f;
	color[1] = 6.0f * v2 * v2 * (1.0f - v2);
	color[2] = 5.5f * v3 * (1.0f - v3) * (1.0f - v3);

	// alternative: classical jet colormap
	// color[0] = 1.5 - 4.0 * fabs(v - 0.75);
	// color[1] = 1.5 - 4.0 * fabs(v - 0.5 );
	// color[2] = 1.5 - 4.0 * fabs(v - 0.25);
}

// function to fill the data table of the nPoints particles positions, speeds, colors and transparency and the coord table with the nPoints particles positions used to draw;
// data[i][0] == coord[i][0] && data[i][1] == coord[i][1]
void fillData(GLfloat(* data)[8])
{
	float rmax = 100.0 * sqrtf(2.0f);
// 	for (int i = 0; i < NPTS; i++) {
// 		data[i][0] = rand() * 200.0 / RAND_MAX - 100.0; // x (rand between -100 and 100)
// 		data[i][1] = rand() * 200.0 / RAND_MAX - 100.0; // y (rand between -100 and 100)
// 		double r = sqrt(data[i][0] * data[i][0] + data[i][1] * data[i][1]);
// 		data[i][2] = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
// 		data[i][3] = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
// 		colormap(r / rmax, &data[i][4]); // fill color
// 		data[i][7] = 0.8f; // transparency
// 	}
	double x_lim[2] = {-1.0, 1.0};
	for (int i = 0; i < NPTS; i++) {
	  double coord[2] = {data[i][0], data[i][1]};
	  double values[2] = {0.0,0.0};
	  double mass = 1.0;
	  double density = 1.0;
	  
	  init1DSegmentWithParticles(x_lim, coord, values, &mass, &density, NPTS, i, 1);
	  double r = sqrt(data[i][0] * data[i][0] + data[i][1] * data[i][1]);
	  data[i][2] = 0.0; 
	  data[i][3] = 0.0;
	  colormap(r / rmax, &data[i][4]); // fill color
	  data[i][7] = 0.8f; // transparency
	}
}

int main()
{
	GLfloat(*data)[8] = malloc(sizeof(data[0]) * NPTS);
	CHECK_MALLOC(data);
	// Seed the random
	time_t seed = time(NULL);
	//printf(" %u \n", seed);
	srand(seed);
	fillData(data);
	
	// Creation of neighborhoods
	double timestep = 0.0;
	double maxspeed = 0.0;
	neighborhood_options* options = neighborhood_options_init(timestep, maxspeed);
	neighborhood* nh = options->nh;
	neighborhood_update(options, nh, data, 0);
	
	// Creation of set of particles, for which to each of them are associated their coordinates, values of quantity they carry and their neighbourhood
	mySingleParticle* my_array_of_particles = create_array_of_particles(NPTS, 1, nh);
	// Compute the derivatives of the quantity carried by the particles 
	computeDerivatiesAllParticles(my_array_of_particles, NPTS);

	neighborhood_options_delete(options,nh);

// 	double timestep = 0.5;
// 	double maxspeed = 1;
// 	neighborhood_options* options = neighborhood_options_init(timestep, maxspeed);
// 	neighborhood* nh = options->nh;
// 	int number_of_iterations = 10;
// 	for (int iterations = 0; iterations < number_of_iterations;iterations++) {
// 		if(iterations)
// 			bouncyrandomupdate(data, timestep, options->half_length, maxspeed);
// 		neighborhood_update(options, nh, data, iterations);
// 		//kernel(data, nh, kh);
// 	}
// 	neighborhood_options_delete(options,nh);

	free(data);
	free(my_array_of_particles);
	return EXIT_SUCCESS;
}
