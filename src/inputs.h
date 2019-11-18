#include "draw_tools.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

GLfloat random_gauss(GLfloat mu, GLfloat sigma);

void random_uniform_points(GLfloat coord[][2], GLsizei n,
                           GLfloat min[2], GLfloat max[2]);

void random_points(GLfloat coord[][2], GLsizei n);

void random_polygon(GLfloat coord[][2], GLsizei n, int nSmooth);