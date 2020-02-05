#include "BOV.h"

int main()
{
	bov_window_t* window = bov_window_new(800, 800, "Tutorial 1");
	bov_window_set_color(window, (GLfloat[]){0.9f, 0.85f, 0.8f, 1.0f});

	const GLsizei nPoints = 500;
	GLfloat (*data)[8] = malloc(sizeof(data[0])*nPoints);
	for(int i=0; i<nPoints; i++) {
		data[i][0] = rand() * 2.0 / RAND_MAX - 1.0;
		data[i][1] = rand() * 2.0 / RAND_MAX - 1.0;
	}

	bov_points_t *particles = bov_particles_new(data, nPoints, GL_STATIC_DRAW);
	// bov_points_set_width(particles, 0.0002);

	while(!bov_window_should_close(window)){
		bov_particles_draw(window, particles, 0, nPoints);

		bov_window_update(window);
	}

	bov_points_delete(particles);
	free(data);
	bov_window_delete(window);

	return EXIT_SUCCESS;
}