#ifndef PRINT_PARTICULES_H
#define PRINT_PARTICULES_H

#include "BOV.h"
#include <math.h>
#include "particle.h"

typedef struct Animation Animation;

struct Animation {
	bov_window_t* window;
	bov_points_t* particles;
	double timeout;
	int N;
	bov_points_t* grid;
};

Animation* Animation_new(int N, double timeout,Grid* grid);
void Animation_free(Animation* animation);

void fillData(GLfloat(*data)[8], Particle** particles, int N);
bov_points_t * load_Grid(Grid* grid);
void colormap_cell(Particle* p, float color[3]);
void colormap_uni_color(float color[3]);
void colours_neighbors(GLfloat(*data)[8], Particle** particles, int index);

void display_particles(Particle** particles, Animation* animation,bool end);



#endif
