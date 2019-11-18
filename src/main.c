#include "inputs.h"
#include <time.h>


int main()
{
	// give a bit of entropy for the seed of rand()
	// or it will always be the same sequence
	int seed = (int) time(NULL);
	srand(seed);

	// we print the seed so you can get the distribution of points back
	printf("The seed is %d\n", seed);

	window_t* window = window_new(800,800, "Tutorial 1");
	window_set_color(window, (GLfloat[]){0.7f, 0.7f, 0.7f, 1.0f});

	const GLsizei nPoints = 500;
	GLfloat (*coord)[2] = malloc(sizeof(coord[0])*nPoints);
#if 0 // put 1 for random polygon
	random_polygon(coord, nPoints, 4);
#else
	random_points(coord, nPoints);
#endif

	points_t *coordDraw = points_new(coord, nPoints, GL_STATIC_DRAW);
	points_set_color(coordDraw, (GLfloat[4]){0.0, 0.0, 0.0, 1.0});
	points_set_outline_color(coordDraw, (GLfloat[4]){0.3, 0.12, 0.0, 0.1});

	while(!window_should_close(window)){
		points_set_width(coordDraw, 0.005);
		points_set_outline_width(coordDraw, 0.003);
		line_loop_draw(window, coordDraw, 0, nPoints);

		points_set_width(coordDraw, 0.003);
		points_set_outline_width(coordDraw, -1.);
		points_draw(window, coordDraw, 0, nPoints);


		window_update(window);
	}

	points_delete(coordDraw);
	free(coord);
	window_delete(window);

	return EXIT_SUCCESS;
}