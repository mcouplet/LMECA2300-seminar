#ifndef PARTICULES_H
#define PARTICULES_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#include "utils.h"

typedef struct Cell Cell;
typedef struct Grid Grid;
typedef struct Particle Particle;
typedef struct Verlet Verlet;

struct Cell {
	int i,j;
	List* neighboorCells;
	List* particles;
	bool visited; //for improved algorithm
};

struct Grid {
	int nCellx;
	int nCelly;
	Cell** cells;
	double h;    //length of a Cell
	double left;
	double right;
	double top;
	double bottom;
};

struct Particle {
	int index;
	xy* pos;
	xy* v;
	double density;
	double e;

	Cell* cell;
	List* neighborhood;
	List* potential_neighborhood;
};

struct Verlet {
	bool verlet;
	double kh;
	double L;
	int T;
};

Grid* Grid_new(double x1, double x2, double y1, double y2, double kh);
void Grid_free(Grid* grid);
void Cell_free(Cell* cell);

void Particle_set(Particle* p,double x, double y, double vx, double vy, double d, double e,int index);
void Particle_free(Particle* p,int N);

Verlet* Verlet_new(double kh, double L, int T);

//update of the particule's location in cells
Cell* localise_particle(Grid* grid, Particle * p);
void updateCells(Grid* grid, Particle* particles, int N);

//update of the neighborhood
void add_neighbors_from_cell(Particle* p, Cell* cell, double r);
void add_neighbors_from_cells(Grid* grid, Particle* p);
void update_from_potential_neighbors(Particle* particles, int N, double r);
void update_neighborhoods(Grid* grid, Particle* particles, int N, int iter, Verlet* verlet);

//update of the neighborhood improved
void update_neighborhoods_improved(Grid* grid);

//build random particules
Particle* build_Particles(int N,double L);

void Grid_reset(Grid* grid);
void Particles_reset(Particle* particles, int N, int iter, Verlet* verlet);

#endif