#include "Particule.h"


Grid* Grid_new(double x1, double x2, double y1,double y2,double h){
	int nCellx = ceil((x2-x1) / h);
	int nCelly = ceil((y2-y1) / h);
	printf("(%d,%d)\n", nCellx, nCelly);
	//build the grid
	Cell** cells = (Cell**)malloc(nCellx * sizeof(Cell*));
	Cell* tab = (Cell*)malloc(nCellx * nCelly * sizeof(Cell));
	for (int i = 0;i < nCellx;i++)
		cells[i] = &tab[i*nCelly];
	
	
	//build the link between the Cells
	for (int i = 0;i < nCellx;i++){
		for (int j = 0;j < nCelly;j++){
			cells[i][j].i = i;cells[i][j].j = j;
			cells[i][j].neighboorCells = List_new();
			cells[i][j].particles = List_new();
			cells[i][j].visited = false;
			for (int k = -1;k <= 1;k++){
				for (int l = -1;l <= 1;l++){
					if (k != 0 || l != 0)
						if (i + k >= 0 && i + k < nCellx && j + l >= 0 && j + l < nCelly)
							List_append(cells[i][j].neighboorCells,&cells[i + k][j + l]);
				}
			}
		}
	}
	Grid* grid = (Grid*)malloc(sizeof(Grid));
	grid->left = x1;grid->right= x1+ nCellx*h;
	grid->bottom = y1;grid->top = y1 + nCelly * h;
	grid->nCellx = nCellx;
	grid->nCelly = nCelly;
	grid->h = h;
	grid->cells = cells;
	return grid;
}

void Grid_free(Grid* grid)
{
	Cell** cells = grid->cells;
	for (int i = 0;i < grid->nCellx;i++){
		for (int j = 0;j < grid->nCelly;j++){
			Cell_free(&cells[i][j]);
		}
	}
	free(cells[0]);
	free(cells);
	free(grid);
}

void Cell_free(Cell* cell)
{
	List_free(cell->neighboorCells, NULL);
	List_free(cell->particles, NULL);
}

void Particle_set(Particle* p,double x,double y,double vx,double vy,double d,double e,int index)
{
	p->index = index;
	p->pos = xy_new(x,y);
	p->v = xy_new(vx, vy);
	p->density = d;
	p->e = e;
	p->cell = NULL;
	p->neighborhood = List_new();
	p->potential_neighborhood = List_new();
}
void Particle_free(Particle* particles,int N)
{
	for (int i = 0;i < N;i++){
		free(particles[i].pos);
		free(particles[i].v);
		List_free(particles[i].neighborhood, NULL);
		List_free(particles[i].potential_neighborhood, NULL);
	}
	free(particles);
}

Verlet* Verlet_new(double kh,double L,int T)
{
	Verlet* v = (Verlet*)malloc(sizeof(Verlet));
	v->kh = kh;
	v->L = L;
	v->T = T;
	return v;
}
///////////////////////////update cells///////////////////////////
Cell* localise_particle(Grid* grid, Particle * p)
{
	int i = floor((p->pos->x - grid->left)/grid->h);
	int j = floor((p->pos->y - grid->bottom) / grid->h);
	return &grid->cells[i][j];
}

void updateCells(Grid* grid, Particle* particles,int N)
{
	//clean the grid before update
	Grid_reset(grid);
	for (int i = 0;i < N;i++){
		Cell* cell = localise_particle(grid,&particles[i]);
		particles[i].cell = cell;
		List_append(cell->particles, &particles[i]);
	}
}
/////////////////////////update neighborhood////////////////////////
void add_neighbors_from_cell(Particle* p, Cell* cell , double r)
{
	ListNode *node = cell->particles->head;
	while (node != NULL) {
		Particle* q = (Particle*)node->v;
		if ((p->index != q->index) && pow(p->pos->x - q->pos->x, 2) + pow(p->pos->y - q->pos->y, 2) <= pow(r,2))
				List_append(p->neighborhood, q);
		node = node->next;
	}
}
void add_neighbors_from_cells(Grid* grid, Particle* p)
{
	add_neighbors_from_cell(p, p->cell, grid->h);
	ListNode *node = p->cell->neighboorCells->head;
	while (node != NULL) {
		Cell* cell = (Cell*)node->v;
		add_neighbors_from_cell(p, cell, grid->h);
		node = node->next;
	}
}

//take among the potentail neighbors, the real neighbors
void update_from_potential_neighbors(Particle* particles, int N, double r)
{
	for (int i = 0;i < N;i++)
	{
		Particle* p = &particles[i];
		ListNode *node = p->potential_neighborhood->head;
		while (node != NULL){
			Particle* q = (Particle*)node->v;
			if (pow(p->pos->x - q->pos->x, 2) + pow(p->pos->y - q->pos->y, 2) <= pow(r, 2))
				List_append(p->neighborhood, q);
			node = node->next;
		}
	}
}
void update_neighborhoods(Grid* grid, Particle* particles, int N,int iter, Verlet* verlet)
{
	//clean the particles before update
	Particles_reset(particles,N,iter,verlet);
	if (verlet==NULL)
	{
		//update_neighborhoods_improved(grid);
		for (int i = 0;i < N;i++)
			add_neighbors_from_cells(grid, &particles[i]);
	}
	else
	{
		if (iter%verlet->T == 0) {
			//update_neighborhoods_improved(grid);
			for (int i = 0;i < N;i++) {
				add_neighbors_from_cells(grid, &particles[i]);
				List* l = particles[i].potential_neighborhood;
				particles[i].potential_neighborhood = particles[i].neighborhood;
				particles[i].neighborhood = l;
			}
		}
		update_from_potential_neighbors(particles, N, verlet->kh);
	}
}

//////////////////update neighborhood-improved algorithm///////////////

//compute the neighbors inside a cell
void update_neighbors_1(Cell* cell,double r)
{
	ListNode *node = cell->particles->head;
	while (node != NULL) {
		Particle* p = (Particle*)node->v;

		ListNode *node2 = node->next;
		while (node2 != NULL)
		{
			Particle* q = (Particle*)node2->v;
			if (pow(p->pos->x - q->pos->x, 2) + pow(p->pos->y - q->pos->y, 2) <= pow(r, 2))
			{
				List_append(p->neighborhood, q);
				List_append(q->neighborhood, p);
			}
			node2 = node2->next;
		}

		node = node->next;
	}
}

//compute the neighbouring relation between 2 cells
void neighbors_in_2_cells(Cell* cell1, Cell* cell2, double r)
{
	ListNode *node1 = cell1->particles->head;
	while (node1 != NULL) {
		Particle* p1 = (Particle*)node1->v;
		
		ListNode *node2 = cell2->particles->head;
		while (node2 != NULL) {
			Particle* p2 = (Particle*)node2->v;
			if (pow(p1->pos->x - p2->pos->x, 2) + pow(p1->pos->y - p2->pos->y, 2) <= pow(r, 2))
			{
				List_append(p1->neighborhood, p2);
				List_append(p2->neighborhood, p1);
			}
			node2 = node2->next;
		}

		node1 = node1->next;
	}
}
//compute the neighbors of particules of cell with the neighboring cells
void update_neighbors_2(Cell* cell, double r)
{
	ListNode *node = cell->neighboorCells->head;
	while (node != NULL) {
		Cell* neighborCell = (Cell*)node->v;
		if (!neighborCell->visited)
		{
			neighbors_in_2_cells(cell, neighborCell, r);
		}
		node = node->next;
	}
}
void update_neighborhoods_improved(Grid* grid)
{
	Cell** cells = grid->cells;
	for (int i = 0;i < grid->nCellx;i++)
	{
		for (int j = 0;j < grid->nCelly;j++)
		{
			update_neighbors_1(&cells[i][j], grid->h);
			update_neighbors_2(&cells[i][j], grid->h);
			cells[i][j].visited = true;
		}
	}
}


////////////////////////////////////////////////////

//generate random particles
Particle* build_Particles(int N,double L)
{
	Particle* particles = (Particle*) malloc(N*sizeof(Particle));
	for (int i = 0; i < N; i++) {
		double x = rand_interval(-L, L);
		double y = rand_interval(-L, L);
		Particle_set(&particles[i],x, y, 0, 0, 0, 0, i);
	}
	return particles;
}

//reset the list of particles inside each Cell
void Grid_reset(Grid* grid)
{
	for (int i = 0;i < grid->nCellx;i++)
	{
		for (int j = 0;j < grid->nCelly;j++)
		{
			Cell* cell = &grid->cells[i][j];
			List_free(cell->particles, NULL);
			
			cell->particles = List_new();
			cell->visited = false;
		}
	}
}
//reset neighborhood of each particles
void Particles_reset(Particle* particles,int N,int iter,Verlet* verlet)
{
	for (int i = 0;i < N;i++)
	{
		Particle* p = &particles[i];
		List_free(p->neighborhood,NULL);
		p->neighborhood = List_new();
		if (verlet != NULL && iter%verlet->T==0) {
			List_free(p->potential_neighborhood, NULL);
			p->potential_neighborhood = List_new();
		}
	}
}


