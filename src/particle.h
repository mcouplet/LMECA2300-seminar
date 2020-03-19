#ifndef PARTICLE_H
#define PARTICLE_H

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
	List* neighboring_cells;
	List* particles;
	bool visited; // for improved algorithm
};

struct Grid {
	int nCellx;
	int nCelly;
	Cell*** cells;	// 2D array of cell pointers
	double h; // side length of a cell
	double left;	// x-coordinate of left side
	double right;	// x-coordinate of right side
	double top;		// y-coordinate of top side
	double bottom;	// y-coordinate of bottom side
};

struct Particle {
	int index;
	xy* pos; // position
	xy* v; // velocity
	double density;
	double e; // ???

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

void Cell_free(Cell* cell); // Cell destructor

Grid* Grid_new(double x1, double x2, double y1, double y2, double kh); // Grid constructor
void Grid_free(Grid* grid); // Grid destructor

//Particle* Particle_new(int index, xy* pos, xy* v, double density, double e)
void Particle_set(Particle* p, double x, double y, double vx, double vy, double d, double e, int index); // Particle setter

void free_particles(Particle** p, int N); // destroy array of particles

Verlet* Verlet_new(double kh, double L, int T);

// Update of the particles' locations in cells
Cell* localize_particle(Grid* grid, Particle * p);
void update_cells(Grid* grid, Particle** particles, int N);

// Update of the neighborhood
void add_neighbors_from_cell(Particle* p, Cell* cell, double r);
void add_neighbors_from_cells(Grid* grid, Particle* p);
void update_from_potential_neighbors(Particle** particles, int N, double r);
void update_neighborhoods(Grid* grid, Particle** particles, int N, int iter, Verlet* verlet);

// Improved update of the neighborhood
void update_neighborhoods_improved(Grid* grid);

// Build random particles
Particle** build_particles(int N, double L);

void reset_grid(Grid* grid);
void reset_particles(Particle** particles, int N, int iter, Verlet* verlet);

#endif
