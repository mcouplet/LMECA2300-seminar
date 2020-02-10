#include "BOV.h"
#include <math.h>


void colormap(float v, float color[3])
{
	float v1 = 3.5*(v-0.7);
	float v2 = 1.25*v;
	float v3 = fminf(0.5,v)*2.0;
    
    color[0] = 0.8f*(-v1*v1+1.0f);
    color[1] = 0.8f*(6.0f*v2*v2*(1.0f-v2));
    color[2] = 0.8f*(5.5f*v3*(1.0f-v3)*(1.0f-v3));
}


int main()
{
	bov_window_t* window = bov_window_new(800, 800, "Tutorial 1");
	bov_window_set_color(window, (GLfloat[]){0.9f, 0.85f, 0.8f, 1.0f});

	const GLsizei nPoints = 10000;
	GLfloat (*data)[8] = malloc(sizeof(data[0])*nPoints);

	float rmax = sqrtf(2.0f);
	for(int i=0; i<nPoints; i++) {
		data[i][0] = rand() * 2.0 / RAND_MAX - 1.0; // x
		data[i][1] = rand() * 2.0 / RAND_MAX - 1.0; // y
		data[i][2] = 0; // speed x (not used by default visualization)
		data[i][3] = 0; // speed y (not used by default visualization)
		float r = sqrt(data[i][0]*data[i][0] + data[i][1]*data[i][1]);
		colormap(r/rmax, &data[i][4]); // fill color
		data[i][7] = 0.8f; // transparency

	}

	bov_points_t *particles = bov_particles_new(data, nPoints, GL_STATIC_DRAW);
	bov_points_set_width(particles, 0.01);
	bov_points_set_outline_width(particles, 0.0025);

	while(!bov_window_should_close(window)){
		bov_particles_draw(window, particles);
		// bov_points_draw(window, particles, 0, BOV_TILL_END);

		bov_window_update(window);
	}

	bov_points_delete(particles);
	free(data);
	bov_window_delete(window);

	return EXIT_SUCCESS;
}