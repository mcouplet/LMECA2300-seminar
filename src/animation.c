#include "draw_tools.h"
#include <math.h>

const float linewidth = 0.03;
const float transition_time = 1.0;

// transition from a point p to another point
void transition(points_t* pointset, float coord[20], GLuint p, float x)
{
	// we don't want p to move with the scaling
	points_set_pos(pointset, (GLfloat[2]){coord[2*p]*(1-x), coord[2*p+1]*(1-x)});

	// scale everything
	points_scale(pointset, (GLfloat[2]){x, x});
}

// reseting the parameters of the object
void reset_params(points_t* pointset)
{
	points_set_color(pointset, (float[4]){1.0, 0.6, 0.3, 1.0});
	points_set_pos(pointset, (GLfloat[2]){0, 0});
	points_set_width(pointset, linewidth);
	points_scale(pointset, (GLfloat[2]){1.0, 1.0});
}


int main(int argc, char *argv[])
{
	window_t* window = window_new(0,0, argv[0]);

	float coord[20] = {-0.2, -0.4,
                       -0.5, 0.0,
                        0.3, 0,
                        0.0, -0.9,
                       -0.8, -0.6,
                       -1.0, 0,
                       -0.7, 0.6,
                        0.0, 0.8,
                        0.7, 0.6,
                        0.5, -0.6};

	points_t* pointset = points_new(coord, 10, GL_STATIC_DRAW);

	// by passing a null pointer, we only reserve space...
	order_t* full_order = order_new(NULL, 10, GL_DYNAMIC_DRAW);
	order_t* transition_order = order_new(NULL, 2, GL_DYNAMIC_DRAW);

	// you can change the order in indices for a less obvious one...
	GLuint indices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0};
	int len = sizeof(indices)/sizeof(indices[0]);

	// a grey background
	window_set_color(window, (GLfloat[4]){0.3, 0.3, 0.3, 1});

	for (int i=0; i<len-1; i++) {
		double tbegin = window_get_time(window);
		double tnow = tbegin;

		GLuint transition_pair[2] = {indices[i], indices[i+1]};


		transition_order = order_update(transition_order, transition_pair, 2, GL_DYNAMIC_DRAW);

		while(tnow - tbegin < transition_time) {
			if(window_should_close(window))
				goto end_of_program; // break all the loop (only valid use of goto)

			reset_params(pointset);
			line_strip_draw_with_order(window, pointset, full_order);
			transition(pointset, coord, indices[i], (tnow-tbegin)/transition_time);
			lines_draw_with_order(window, pointset, transition_order);

			window_update(window);
			tnow = window_get_time(window);
		}

		full_order = order_update(full_order, indices, i+2, GL_DYNAMIC_DRAW);
	}

	// we want to keep the window open with everything displayed...
	reset_params(pointset);
	while(!window_should_close(window)) {
		line_strip_draw_with_order(window, pointset, full_order);
		window_update_and_wait_events(window);
	}

end_of_program:

	printf("Ended correctly - %.2f second\n", window->wtime);

	points_delete(pointset);
	order_delete(full_order);
	order_delete(transition_order);
	window_delete(window);

	return EXIT_SUCCESS;
}
