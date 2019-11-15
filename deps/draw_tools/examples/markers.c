 /*************************************************************************
  * Animation example program using Draw_tools, a wrapper around OpenGL and
  * GLFW (www.glfw.org) to draw simple 2D graphics.
  *------------------------------------------------------------------------
  * Copyright (c) 2019-2020 Célestin Marot <marotcelestin@gmail.com>
  *
  * This software is provided 'as-is', without any express or implied
  * warranty. In no event will the authors be held liable for any damages
  * arising from the use of this software.
  *
  * Permission is granted to anyone to use this software for any purpose,
  * including commercial applications, and to alter it and redistribute it
  * freely, subject to the following restrictions:
  *
  * 1. The origin of this software must not be misrepresented; you must not
  *    claim that you wrote the original software. If you use this software
  *    in a product, an acknowledgment in the product documentation would
  *    be appreciated but is not required.
  *
  * 2. Altered source versions must be plainly marked as such, and must not
  *    be misrepresented as being the original software.
  *
  * 3. This notice may not be removed or altered from any source
  *    distribution.
  *
  *************************************************************************/

#include "draw_tools.h"
#include <math.h>

// Maybe there will be more marker types available in the future.
// At the moment, markers value wraps around at 25
// We cannot have more than 100 markers without chaning the way the
// string is created (putting %7.3f and lengthen the string)
#define NUMMARKERS 25


int main(int argc, char *argv[])
{
	/* Actually, this is not an example at all, you shouldn't draw the
	 * same f**king point 250 time, each time with a different draw
	 * call.
	 *
	 * This is more of a tool. When you launch this program, you will
	 * be able to see which marker value correspond to which shape
	 * and choose the shape that fits your need accordingly.
	 *
	 * REMINDER: DO NOT MAKE THAT MANY DRAW CALL IN YOUR CODE !!!!
	 */

	window_t* window = window_new(0,0, argv[0]);
	window_set_color(window, (float[4]){1.0, 0.8, 0.5, 1});

	const float pointWidth = 1.0/NUMMARKERS;

	points_t* single_point = points_new((float[2]){0,0}, 1, GL_STATIC_DRAW);
	points_set_outline_color(single_point, (float[4]){0.3, 0.3, 0.3, 1});
	points_set_outline_width(single_point, pointWidth*0.2);
	points_set_width(single_point, pointWidth);

	text_t* marker_text;
	{
		/* first, we will generate all the text values */
		// we will have 6 charater per number (2digit, one dot, 3 digit)
		// then 2 space, for each marker.
		// then we have 10 \n at the end of each line
		char string[9*(NUMMARKERS*8 + 9)+1];
		int cur = 0;

		for(int j=0; j<9; j++) {
			for(int i=0; i<NUMMARKERS; i++) {
				sprintf(string+cur, " %6.3f ", i+j*0.12493);
				cur += 8;
			}

			// replace the 2 trailing space with 4 '\n'
			sprintf(string+cur-1, "\n\n\n\n\n\n\n\n\n\n");
			cur+= 9;
		}


		marker_text = text_new((unsigned char *) string,
			                           GL_STATIC_DRAW);
		text_set_fontsize(marker_text, pointWidth*0.5);
		text_set_pos(marker_text, (GLfloat[2]){-1, 1-pointWidth});
	}
	

	while(!window_should_close(window)){
		double wtime = window_get_time(window);
		double fract02 = modf(0.2*wtime, &wtime);
		wtime = fabs(2*fract02-1.0);

		// we change the color over time
		points_set_color(single_point, (float[4]) {
		                 sin(0.11*wtime)*0.5+0.5,
		                 sin(0.7*wtime)*0.5+0.5,
		                 sin(0.67*wtime)*0.5+0.5,
		                 1});

		for(int i=0; i<NUMMARKERS; i++) {

			GLfloat pos[2] = {pointWidth-1.0+2.0*pointWidth*i,
			                 1.0-2.5*pointWidth};

			for(int j=0; j<9; j++) {
				points_set_pos(single_point, pos);
				points_set_marker(single_point, i+j*0.12493);
				points_draw(window, single_point, 0, 1);
				
				pos[1] -= 5*pointWidth;
			}

			points_set_pos(single_point, pos);
			points_set_marker(single_point, i+wtime);
			points_draw(window, single_point, 0, 1);
		}

		text_draw(window, marker_text);

		window_update(window);
	}

	printf("Ended correctly");

	points_delete(single_point);
	text_delete(marker_text);
	window_delete(window);

	return EXIT_SUCCESS;
}
