
#include "checkDerivatives.h"

singleParticleDerivatives* initialize_particle_derivatives(int size_values)
{
  singleParticleDerivatives* my_deriv = malloc(sizeof(singleParticleDerivatives));
  
  if (size_values == 1) { // scalar quantity
      my_deriv->divergence = NULL; //  divergence of a scalar has no meaning (= gradient)
      my_deriv->gradient = (double*) calloc(2,sizeof(double));
      my_deriv->laplacian = (double*) calloc(1, sizeof(double));
  }
  else if (size_values == 2) { // vector quantity
      my_deriv->divergence = (double*) calloc(1, sizeof(double));
      my_deriv->gradient = (double*) calloc(4,sizeof(double));
      my_deriv->laplacian = (double*) calloc(2, sizeof(double));
  }
  else {
    printf("---------------  The size of the quantity for which you try to compute the derivatives is not correct ------------------------");
    return 0;
  }
  
  return my_deriv;
}

mySingleParticle* create_array_of_particles(int nbParticles, int size_values, neighborhood* nh) 
{
  
  mySingleParticle* my_particles = malloc(nbParticles*sizeof(mySingleParticle));
  
//   double x_lim[2] = {-1.0, 1.0};
//   for (int i = 0; i < nbParticles; i++) {
//     my_particles[i].coordinates = calloc(2,sizeof(double));
//     my_particles[i].values = calloc(size_values,sizeof(double));
//     my_particles[i].mass = 0.0;
//     my_particles[i].density = 0.0;
//     my_particles[i].size_values = size_values;
//     my_particles[i].particle_neighbours = neighborhood_options_init(0.0, 0.0);
//     my_particles[i].particle_neighbours->nh = &(nh[i]);
//     my_particles[i].particle_derivatives = initialize_particle_derivatives(size_values);
//     
//     init1DSegmentWithParticles(x_lim, my_particles[i].coordinates, my_particles[i].values, &(my_particles[i].mass), &(my_particles[i].density), nbParticles, i, size_values);  
//   }
  double x_lim[4] = {-1.0, 1.0, -1.0, 1.0};
  int nbPart_x = sqrt(nbParticles);
  int nbPart_y = nbPart_x;
  int i = 0;
  for (int ind_y = 0; ind_y < nbPart_y; ind_y++) {
      for (int ind_x = 0; ind_x < nbPart_x; ind_x++) {
	my_particles[i].coordinates = calloc(2,sizeof(double));
	my_particles[i].values = calloc(size_values,sizeof(double));
	my_particles[i].mass = 0.0;
	my_particles[i].density = 0.0;
	my_particles[i].size_values = size_values;
	my_particles[i].particle_neighbours = neighborhood_options_init(0.0, 0.0);
	my_particles[i].particle_neighbours->nh = &(nh[i]);
	my_particles[i].particle_derivatives = initialize_particle_derivatives(size_values);
    
	initSquareWithParticles(x_lim, my_particles[i].coordinates, my_particles[i].values, &(my_particles[i].mass), &(my_particles[i].density), nbPart_x, ind_x, ind_y, size_values);  
	i++;
      }
   }
  return my_particles;
}

void initSquareWithParticles(double* x_lim, double* coord, double* values, double* mass, double* density, int nb_particles, int index_x, int index_y, int size_values)// double (*myFun)(double*)) 
{
  double x_min = x_lim[0];
  double x_max = x_lim[1];
  double y_min = x_lim[2];
  double y_max = x_lim[3];
  if (x_min > x_max) {
      double x_temp = x_min;
      x_min = x_max;
      x_max = x_temp;
  }
  if (y_min > y_max) {
      double y_temp = y_min;
      y_min = y_max;
      y_max = y_temp;
  }
  double length_x = fabs(x_max-x_min);
  double delta_x = length_x / ((double)nb_particles - 1.0);
  double length_y = fabs(y_max-y_min);
  double delta_y = length_y / ((double)nb_particles - 1.0);

  coord[0] = x_min + index_x*delta_x;
  coord[1] = y_min + index_y*delta_y;

  // function values
  values[0] = myFunctionToDerive(coord);//(myFun)(&x_coord);
  if (size_values > 1) values[0] = 0.0; // 1-D segment so don't care about the second dimension
  *mass = 1.0;
  *density = *mass / (delta_x*delta_y);//WARNING: only valid for equispaced grid of particle! (double)nb_particles / length_x;//1.0; 
}

void init1DSegmentWithParticles(double* x_lim, double* coord, double* values, double* mass, double* density, int nb_particles, int index_part, int size_values)// double (*myFun)(double*)) 
{
  double x_min = x_lim[0];
  double x_max = x_lim[1];
  if (x_min > x_max) {
      double x_temp = x_min;
      x_min = x_max;
      x_max = x_temp;
  }
  double length = fabs(x_max-x_min);
  double delta_x = length / ((double)nb_particles - 1.0);
  double x_coord = x_min + index_part*delta_x; // uniformly distributed points on the segment
  // coordinates
  coord[0] = x_coord ; // X-dim
  coord[1] = 0.0; // 1-D segment so don't care about the second dimension
  // function values
  values[0] = myFunctionToDerive(coord);//(myFun)(&x_coord);
  if (size_values > 1) values[1] = 0.0; // 1-D segment so don't care about the second dimension
  *mass = 1.0;
  *density = (double)nb_particles / length;//1.0; 
}

double myFunctionToDerive(double* x) 
{
   return x[0]*x[0]; 
}

void computeDerivativesOfParticleQuantity(mySingleParticle* myPart, int index_part) {
   
  mySingleParticle* local_part = &(myPart[index_part]);
  double * divergence_part = local_part->particle_derivatives->divergence;
  double * grad_part = local_part->particle_derivatives->gradient;
  double * laplacian_part = local_part->particle_derivatives->laplacian;

  int nNeigh = local_part->particle_neighbours->nh->nNeighbours;
  neighbours* List = local_part->particle_neighbours->nh->list;
  double kh = local_part->particle_neighbours->kh;
  for (int j = 0; j < nNeigh; j++) {
      int index_j = List->index;
      double d_ij = List->distance; // WARNING: distance = 0 if computed this way, why???
      double d_x_ij =  myPart[index_j].coordinates[0] - local_part->coordinates[0];//  data[index_node2][0] - data[i][0];
//       double d_ij = fabs(d_x_ij); // only valid in 1-D, should use List->distance instead but = 0 for the moment (see line above)
      double d_y_ij = myPart[index_j].coordinates[1] - local_part->coordinates[1];//data[index_node2][1] - data[i][1];
      
      /*
	You can choose here the desired kernel function for your code.
	*/
      
      double weight_x = grad_w_cubic(d_ij, kh, d_x_ij);
      double weight_y = grad_w_cubic(d_ij, kh, d_y_ij);

//       double weight_x = grad_w_lucy(d_ij, kh, d_x_ij);
//       double weight_y = grad_w_lucy(d_ij, kh, d_y_ij);
      
//       double weight_x = grad_w_newquartic(d_ij, kh, d_x_ij);
//       double weight_y = grad_w_newquartic(d_ij, kh, d_y_ij);
      
//       double weight_x = grad_w_quinticspline(d_ij, kh, d_x_ij);
//       double weight_y = grad_w_quinticspline(d_ij, kh, d_y_ij);
      
      // --- Divergence ---
      if (myPart->size_values == 2) divergence_part[0] += (1.0 / local_part->density) * ( (myPart[index_j].values[0] - local_part->values[0]) * weight_x 
							  + (myPart[index_j].values[1] - local_part->values[1]) * weight_y ) * myPart[index_j].mass;
      // --- Gradient & Laplacian ---
      for (int k=0; k < myPart->size_values; k++) {
	// d/dx value_k
	grad_part[2*k] += -local_part->density * myPart[index_j].mass * ((local_part->values[k] / (local_part->density*local_part->density)) 
								      + (myPart[index_j].values[k] / (myPart[index_j].density*myPart[index_j].density))) * weight_x;							      
	// d/dy value_k
	grad_part[2*k+1] += -local_part->density * myPart[index_j].mass * ((local_part->values[k] / (local_part->density*local_part->density)) 
								      + (myPart[index_j].values[k] / (myPart[index_j].density*myPart[index_j].density))) * weight_y;
	// d^2/dx^2 value_k + d^2/dy^2 value_k
	laplacian_part[k] += 2.0 * (myPart[index_j].mass / myPart[index_j].density) * (local_part->values[k] - myPart[index_j].values[k]) * (d_x_ij * weight_x + d_y_ij * weight_y) /(d_ij*d_ij);
      }
      
      List = List->next;
  }
// Don't free otherwize you loose the pointers to local_part->particle_derivatives
//   free(divergence_part);
//   free(grad_part);
//   free(laplacian_part);

}


void computeDerivatiesAllParticles(mySingleParticle* myPart, int nbParticles) {
  // WARNING: this routine only valid for scalar quantity for the moment
  printf("\ni  (X,Y)			f(x)		Grad_x			Grad_y			Laplacian\n\n");
  for (int i=0; i<nbParticles; i++) {
    
      computeDerivativesOfParticleQuantity(myPart, i);
      double* x = myPart[i].coordinates;
      double* grad = myPart[i].particle_derivatives->gradient;
      double* lapl = myPart[i].particle_derivatives->laplacian;
      double value = myPart[i].values[0];
      if (x[1] > 0.0 && x[1] < 0.1) {
	printf("%d  (%2.3f,%2.3f)		%2.6f        %2.6f		%2.6f		%2.6f\n",i,x[0],x[1],value,grad[0],grad[1],lapl[0]);
      }
	
     }
  
  
  
}  
  
  