#include "draw_tools.h"
#include <math.h>


int main(int argc, char *argv[])
{
	window_t* window = window_new(0, 0, argv[0]);
	window_set_color(window, (GLfloat[4]) {0.8,0.8,0.8,1.0});

	text_t* outline = text_new((unsigned char[]) {"varying outline width"},
		                        GL_STATIC_DRAW);
	text_param_t parameters = {.outlineColor={1,0,0,1},
	                           .pos={-1.0, 0.33},
	                           .fillColor={0},// completely transparent
	                           .height=0.25,
	                           .boldness=0.25,
	                           .outlineWidth=0.5};
	text_set_param(outline, parameters);

	text_t* width = text_new((unsigned char[]) {"varying width"},
		                     GL_STATIC_DRAW);

	text_set_param(width, parameters);
	text_set_outline_width(width, -1.0);
	text_set_color(width, (GLfloat[4]){0.2, 0.2, 0.2, 1});
	text_set_pos(width, (GLfloat[2]){-1.0, 0.0});

	text_t* shift = text_new((unsigned char[]) {"varying outline shift\n(fake 3D)"},
		                     GL_STATIC_DRAW);
	text_set_param(shift, parameters);
	text_set_pos(shift, (GLfloat[2]){-1.0, -0.33});

	text_t* pixel = text_new((unsigned char[]) {"This text is unmoovable and unzoomable."
	                         "\tIts position and its size (height) must be given in pixels"},
	                         GL_STATIC_DRAW);
	text_set_space_type(pixel, PIXEL_SPACE);

	text_t* unzoomable = text_new((unsigned char[]) {".you can't zoom on this point"},
	                              GL_STATIC_DRAW);
	GLfloat pixel64 = 2.0/window_get_yres(window)*64.0; // 64 pixels in screen coordinates
	text_set_pos(unzoomable, (GLfloat[2]){-1.0, 1-1.1*pixel64});
	text_set_height(unzoomable, pixel64); 
	text_set_space_type(unzoomable, UNZOOMABLE_SPACE);

	

	while(!window_should_close(window)){
		double wtime = window_get_time(window);

		text_set_outline_width(outline, 0.7*sin(wtime)+0.5);
		text_draw(window, outline);

		text_set_boldness(width, 0.5*sin(wtime)-0.1);
		text_draw(window, width);

		text_set_outline_shift(shift, (GLfloat[2]){1.0*sin(3*wtime), 1.0*cos(3*wtime)});
		text_draw(window, shift);

		text_draw(window, pixel);
		text_draw(window, unzoomable);

		window_update(window);
	}

	printf("Ended correctly\n");

	text_delete(outline);
	text_delete(width);
	text_delete(shift);

	window_delete(window);

	return EXIT_SUCCESS;
}
