#include "draw_tools.h"
#include <math.h>


int main(int argc, char *argv[])
{
	window_t* window = window_new(1080,640, argv[0]);

	text_t* hello_world = text_new((unsigned char []){"Hello World !"},
		                           GL_STATIC_DRAW);

	while(!window_is_closed(window)){
		text_draw(window, hello_world);

		window_update(window);
	}

	text_delete(hello_world);
	window_delete(window);

	return EXIT_SUCCESS;
}
