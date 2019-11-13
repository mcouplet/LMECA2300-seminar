#include "draw_tools.h"
#include <math.h>

const float transition_time = 1.0;

// transition from a point p to another point
void transition(points_t* diag, float a[2], float b[2], float x)
{
	GLfloat v[2] = {b[0] - a[0], b[1] - a[1]};
	// scale everything
	points_scale(diag, (GLfloat[2]){v[0]*x, v[1]*x});
	points_set_pos(diag, a);
}


int main(int argc, char *argv[])
{
	window_t* window = window_new(0,0, argv[0]);
	
	// a grey background
	window_set_color(window, (GLfloat[4]){0.3, 0.3, 0.3, 1});

	// we define a style for the lines
	points_param_t lineParams = {.fillColor={1.0, 0.6, 0.3, 1.0},
	                             .scale={1.0, 1.0},
	                             .width=0.03};

	GLfloat coord[10][2] = {{-0.2, -0.4},
	                        {-0.5,  0.0},
	                        { 0.3,  0.0},
	                        { 0.0, -0.9},
	                        {-0.8, -0.6},
	                        {-1.0,  0.0},
	                        {-0.7,  0.6},
	                        { 0.0,  0.8},
	                        { 0.7,  0.6},
	                        { 0.5, -0.6}};

	points_t* pointset = points_new((GLfloat*) coord, 10, GL_STATIC_DRAW);
	points_t* diag = points_new((GLfloat[4]){0.0, 0.0, 1.0, 1.0}, 2, GL_STATIC_DRAW);
	points_set_param(pointset, lineParams);
	points_set_param(diag, lineParams);

	for (int i=0; i<10; i++) {
		double tbegin = window_get_time(window);
		double tnow = tbegin;	

		while(tnow - tbegin < transition_time) {
			if(window_should_close(window))
				goto end_of_program; // break all the loop (only valid use of goto)

			line_strip_draw(window, pointset, 0, i+1);

			transition(diag, coord[i], coord[(i+1)%10], (tnow-tbegin)/transition_time);
			lines_draw(window, diag, 0, 2);

			window_update(window);
			tnow = window_get_time(window);
		}
	}

	// we want to keep the window open with everything displayed...
	while(!window_should_close(window)) {
		line_loop_draw(window, pointset, 0, 10);
		window_update_and_wait_events(window);
	}

end_of_program:

	printf("Ended correctly - %.2f second\n", window->wtime);

	points_delete(pointset);
	points_delete(diag);
	window_delete(window);

	return EXIT_SUCCESS;
}
