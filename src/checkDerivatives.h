#ifndef __CHECKDERIVATIES__
#define __CHECKDERIVATIES__

#include "BOV.h"
#include "kernel.h"
#include "neighborhood_search.h"
#include <time.h>
#include <math.h>

// enum fieldNames {
//   Density = 1,
//   Velocity = 2,
//   Pressure = 3,
//   Temperature = 4
// };

typedef struct singleParticleDerivatives {
  double* divergence;
  double* gradient;
  double* laplacian;
} singleParticleDerivatives;

typedef struct mySingleParticle {
  double* coordinates; // x-y
  double* values; // density or velocity (vector!) or pressure or temperature or anything else
  int size_values;
  double mass;
  double density;
  neighborhood_options* particle_neighbours; // all the neighbouring particles infos related to the particle of interest
  singleParticleDerivatives* particle_derivatives; // divergence, gradient and laplacian (depending if the quantity carried by the particle is a scalar or vector)
} mySingleParticle;

// Create an array of particles "mySingleParticle" by initializing the different 
mySingleParticle* create_array_of_particles(int nbParticles, int size_values, neighborhood* nh);

// Initialize the elements of a singleParticleDerivatives
singleParticleDerivatives* initialize_particle_derivatives(int size_values);

// Assign coordinates, values, mass and density to one particle located on a 1-D segment whose limits are in x_lim, with index "index_art". 
void init1DSegmentWithParticles(double* x_lim, double* coord, double* values, double* mass, double* density, int nb_particles, int index_part, int size_values);
void initSquareWithParticles(double* x_lim, double* coord, double* values, double* mass, double* density, int nb_particles, int index_x, int index_y, int size_values);

// Function assigning the values of the quantity carried by each particle based on the coordinates
double myFunctionToDerive(double* x);

// Compute the derivatives (divergence, gradients, laplacian) of a single particle with index "index_part"
// ---> Use for that the functions implemented in "neighborhood_search" and "kernel" 
void computeDerivativesOfParticleQuantity(mySingleParticle* myPart, int index_part);

// Loop on all particles and compute the derivatives of the quantity for each of them
void computeDerivatiesAllParticles(mySingleParticle* myPart, int nbParticles);

#endif