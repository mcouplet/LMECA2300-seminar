#ifndef BOUNDARY_H
#define BOUNDARY_H
#include "particle.h"
#include "kernel.h"
#include "derivatives.h"
void reflective_boundary(GLfloat(*data)[8], GLfloat(*coord)[2], double timestep, double xmin, double xmax, double ymin, double ymax, double iter);
void fill_fictitious_particles(GLfloat(*data_fict)[8], GLfloat(*coord_fict)[2], double minTemp, double maxTemp, float(*sol)[5]);
void repulsive_boundary(GLfloat(*data)[8], GLfloat(*data_fict)[8], GLfloat(*coord)[2], GLfloat(*coord_fict)[2] ,double timestep, double xmin, double xmax, double ymin, double ymax, double time);

#endif
