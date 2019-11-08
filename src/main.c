#include "draw_tools.h"
#include <math.h>


int main(int argc, char *argv[])
{
	window_t* window = window_new(0,0, argv[0]);

	text_t* outline = text_new("varying outline width", GL_STATIC_DRAW);
	text_set_height(outline, 0.25);
	text_set_outline_color(outline, (float[4]){1,0,0,1});
	text_set_boldness(outline, 0.25);
	text_set_pos(outline, -1.0, 0.33);

	text_t* width = text_new("varying width", GL_STATIC_DRAW);
	text_set_height(width, 0.25);
	text_set_outline_color(width, (float[4]){1,0,0,1});
	text_set_pos(width, -1.0, 0.00);

	text_t* shift = text_new("varying outline shift\n(fake 3D)", GL_STATIC_DRAW);
	text_set_height(shift, 0.25);
	text_set_color(shift, (float[4]){1,1,1,1});
	text_set_boldness(shift, 0.25);
	text_set_outline_color(shift, (float[4]){0.1,0.1,0.1,1});
	text_set_outline_width(shift, 0.5); // we don't want the outline to become negative with the shift
	text_set_pos(shift, -1.0, -0.33);

	window_set_color(window, (float[4]) {0.8,0.8,0.8,1.0});

	while(!window_should_close(window)){
		double wtime = window_get_time(window);

		double st = 0.5*sin(wtime);

		text_set_outline_width(outline, st+0.5);
		text_draw(window, outline);

		text_set_boldness(width, st);
		text_draw(window, width);

		text_set_outline_shift(shift, 1.0*sin(3*wtime), 1.0*cos(3*wtime));
		text_draw(window, shift);

		window_update(window);
	}

	printf("Ended correctly - %.2f second\n", window_get_time(window));

	text_delete(outline);
	text_delete(width);
	text_delete(shift);

	window_delete(window);

	return EXIT_SUCCESS;
}
