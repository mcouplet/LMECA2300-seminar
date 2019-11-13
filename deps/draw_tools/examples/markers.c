#include "draw_tools.h"
#include <math.h>


int main(int argc, char *argv[])
{
	window_t* window = window_new(0,0, argv[0]);

	points_t* single_point = points_new((float[2]){0,0}, 1, GL_STATIC_DRAW);
	points_set_outline_width(single_point, 0.02);

	points_set_outline_color(single_point, (float[4]){0.3, 0.3, 0.3, 1});
	window_set_color(window, (float[4]){1.0, 0.8, 0.5, 1});

	while(!window_is_closed(window)){
		double wtime = window_get_time(window);

		// we change the color over time
		points_set_color(single_point, (float[4]) { sin(0.11*wtime)*0.5+0.5, sin(0.7*wtime)*0.5+0.5, sin(0.67*wtime)*0.5+0.5 , 1});
		points_set_width(single_point, 0.1);
		
		points_set_pos(single_point, -0.9, 0.5);
		points_set_marker(single_point, 0);
		points_draw(window, single_point);

		points_set_pos(single_point, -0.6, 0.5);
		points_set_marker(single_point, 1);
		points_draw(window, single_point);

		points_set_pos(single_point, -0.3, 0.5);
		points_set_marker(single_point, 2.1);
		points_draw(window, single_point);

		points_set_pos(single_point, 0.0, 0.5);
		points_set_marker(single_point, 3);
		points_draw(window, single_point);

		points_set_pos(single_point, 0.3, 0.5);
		points_set_marker(single_point, 4);
		points_draw(window, single_point);

		points_set_pos(single_point, 0.6, 0.5);
		points_set_marker(single_point, 5);
		points_draw(window, single_point);

		points_set_pos(single_point, 0.9, 0.5);
		points_set_marker(single_point, 0);
		points_draw(window, single_point);

		points_set_pos(single_point, -0.75, 0.25);
		points_set_marker(single_point, 0.5);
		points_draw(window, single_point);

		points_set_pos(single_point, -0.45, 0.25);
		points_set_marker(single_point, 1.5);
		points_draw(window, single_point);

		points_set_pos(single_point, -0.15, 0.25);
		points_set_marker(single_point, 2.5);
		points_draw(window, single_point);

		points_set_pos(single_point, 0.15, 0.25);
		points_set_marker(single_point, 3.5);
		points_draw(window, single_point);

		points_set_pos(single_point, 0.45, 0.25);
		points_set_marker(single_point, 4.5);
		points_draw(window, single_point);

		points_set_pos(single_point, 0.75, 0.25);
		points_set_marker(single_point, 5.5);
		points_draw(window, single_point);

		points_set_pos(single_point, 0.0, -0.4);
		points_set_marker(single_point, wtime);
		points_set_width(single_point, 0.5);
		points_draw(window, single_point);



		window_update(window);
	}

	printf("Ended correctly - %.2f second\n", window_get_time(window));

	points_delete(single_point);
	window_delete(window);

	return EXIT_SUCCESS;
}
