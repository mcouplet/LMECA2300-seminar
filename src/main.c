#include "draw_tools.h"
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* fonction that outputs a random value, with a
 * probability following a gaussian curve */
static GLfloat random_gauss(GLfloat mu, GLfloat sigma)
{
	GLfloat u1, u2;
	do {
		u1 = (GLfloat)rand() / RAND_MAX;
		u2 = (GLfloat)rand() / RAND_MAX;
	} while (u1 < 1e-6);
	return mu + sigma * sqrtf(-2.0f * logf(u1)) * cosf(2.0f * M_PI * u2);
}


/* fill coord with random coordinates following an uniform distribution */
static void random_uniform_points(GLfloat coord[][2], GLsizei n,
                                  GLfloat min[2], GLfloat max[2])
{
	for (GLsizei i=0; i<n; i++) {
		coord[i][0] = (max[0] - min[0])*rand() / RAND_MAX + min[0];
		coord[i][1] = (max[1] - min[1])*rand() / RAND_MAX + min[1];
	}
}


/* creating random points following a gaussian distribution.
 * around multiple centroid (maximum 6 centroids) which
 * are uniformly*/
static void random_points(GLfloat coord[][2], GLsizei n)
{
	int n_centroids = rand()%6 + 1;
	GLfloat (*centroids)[2] = malloc(sizeof(GLfloat)*n_centroids*2);
	GLfloat (*sigma)[2] = malloc(sizeof(GLfloat)*n_centroids*2);

	GLfloat range = 0.7*(1.0 - 1.0/n_centroids);
	random_uniform_points(centroids, n_centroids,
	                      (GLfloat[2]){-range, -range},
	                      (GLfloat[2]){ range,  range});

	for (GLsizei i=0; i<n_centroids; i++) {
		sigma[i][0] = 0.3f * rand() / RAND_MAX + 0.1f;
		sigma[i][1] = 0.3f * rand() / RAND_MAX + 0.1f;
	}

	for (GLsizei i=0; i<n; i++) {
		for(int j=0; j<2; j++) {
			coord[i][j] = random_gauss(centroids[i%n_centroids][j], sigma[i%n_centroids][j]);
		}
	}

	free(centroids);
	free(sigma);
}


// see https://stackoverflow.com/questions/16542042
static inline GLfloat pseudoangle(GLfloat dx, GLfloat dy) 
{
    GLfloat p = dx/(fabs(dx)+fabs(dy)); // -1 .. 1 increasing with x
    if (dy < 0.0f)
    	return 3.0f + p;  //  2 .. 4 increasing with x
    else
    	return 1.0f - p;      //  0 .. 2 decreasing with x
}


/* compare the angles that two points make 
 *
 * The sorting is not robust, there are a possible floating point
 * errors in the pseudoangle function.
 */
static int compare_angle(const void *a_v, const void *b_v)
{
	GLfloat* a = *(GLfloat(*)[2])a_v;
	GLfloat* b = *(GLfloat(*)[2])b_v;
	
	GLfloat diff = pseudoangle(b[0], b[1]) - pseudoangle(a[0], a[1]);
	return (diff>0) - (diff<0);
}


/* create a random polygon
 * the bigger nSmooth is, the rounder it will be  */
static void random_polygon(GLfloat coord[][2], GLsizei n, int nSmooth)
{
	GLfloat sigmax = (GLfloat)rand() / RAND_MAX;
	GLfloat sigmay = (GLfloat)rand() / RAND_MAX;

	for(GLsizei i=0; i<n; i++) {
		coord[i][0] = random_gauss(0.0f, sigmax);
		coord[i][1] = random_gauss(0.0f, sigmay);
	}

	qsort(coord, n, 2*sizeof(GLfloat), compare_angle);

	// a little bit of smoothing
	for(int smoothing=0; smoothing<nSmooth; smoothing++) {
		// because we do not copy the array and do the smoothing in place
		// we start at a random index :p
		int index = rand()%n;
		for(int i=1; i<n-1; i++) {
			coord[(index+i)%n][0] = (2*coord[(index+i)%n][0] + 
				                       coord[(index+i+n-1)%n][0] +
				                       coord[(index+i+1)%n][0])*0.25f;
			coord[(index+i)%n][1] = (2*coord[(index+i)%n][1] + 
				                       coord[(index+i+n-1)%n][1] +
				                       coord[(index+i+1)%n][1])*0.25f;
		}
	}
}


int main()
{
	// give a bit of entropy for the seed of rand()
	// or it will always be the same sequence
	//srand(time(NULL));

	window_t* window = window_new(800,800, "Tutorial 1");
	window_set_color(window, (GLfloat[]){0.7f, 0.7f, 0.7f, 1.0f});

	const GLsizei nPoints = 100;
	GLfloat (*coord)[2] = malloc(sizeof(coord[0])*nPoints);
	// random_polygon(coord, nPoints, 4);
	random_points(coord, nPoints);

	points_t *coordDraw = points_new(coord, nPoints, GL_STATIC_DRAW);

	while(!window_should_close(window)){
		// points_set_width(coordDraw, 0.001);
		// line_loop_draw(window, coordDraw, 0, nPoints);
		points_set_width(coordDraw, 0.01);
		points_draw(window, coordDraw, 0, nPoints);

		window_update(window);
	}

	points_delete(coordDraw);
	free(coord);
	window_delete(window);

	return EXIT_SUCCESS;
}