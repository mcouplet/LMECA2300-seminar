 /*************************************************************************
  * Points_and_lines example program using BOV, a wrapper around
  * OpenGL and GLFW (www.glfw.org) to draw simple 2D graphics.
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

#include "BOV.h"
#include <math.h>


int main(int argc, char* argv[])
{
	bov_window_t* window = bov_window_new(1080,640, argv[0]);

	bov_text_t* common_text = bov_text_new(
	    (unsigned char[]){"rendering points using"},
	    GL_STATIC_DRAW);
	bov_text_set_pos(common_text, (GLfloat[2]) {-1.0, 0.9});

	bov_text_t* bov_points_draw_label = bov_text_new(
	    (unsigned char[]){
	    "bov_points_draw()"},
	    GL_STATIC_DRAW);

	GLfloat textPos[2] = {-1.0 + 23 * 0.025, 0.9};

	// the default size for a character is 0.025 in width and 0.05 in height
	bov_text_set_pos(bov_points_draw_label, textPos);

	bov_text_t* bov_curve_draw_label = bov_text_new(
	    (unsigned char[]) {"bov_curve_draw()"},
	    GL_STATIC_DRAW);
	bov_text_set_pos(bov_curve_draw_label, textPos);

	bov_text_t* bov_lines_draw_label = bov_text_new(
	    (unsigned char[]) {"bov_lines_draw()"},
	    GL_STATIC_DRAW);
	bov_text_set_pos(bov_lines_draw_label, textPos);

	bov_text_t* bov_lines_draw_with_order_label = bov_text_new(
	    (unsigned char[]) {"bov_lines_draw_with_order()"},
	    GL_STATIC_DRAW);
	bov_text_set_pos(bov_lines_draw_with_order_label, textPos);

	bov_text_t* bov_line_strip_draw_label = bov_text_new(
	    (unsigned char[]) {"bov_line_strip_draw()"},
	    GL_STATIC_DRAW);
	bov_text_set_pos(bov_line_strip_draw_label, textPos);

	bov_text_t* bov_line_strip_draw_with_order_label = bov_text_new(
	    (unsigned char[]) {"bov_line_strip_draw_with_order()"},
	    GL_STATIC_DRAW);
	bov_text_set_pos(bov_line_strip_draw_with_order_label, textPos);

	bov_text_t* bov_line_loop_draw_label = bov_text_new(
	    (unsigned char[]) {"bov_line_loop_draw()"},
	    GL_STATIC_DRAW);
	bov_text_set_pos(bov_line_loop_draw_label, textPos);

	bov_text_t* bov_line_loop_draw_with_order_label = bov_text_new(
	    (unsigned char[]) {"bov_line_loop_draw_with_order()"},
	    GL_STATIC_DRAW);
	bov_text_set_pos(bov_line_loop_draw_with_order_label, textPos);

	bov_points_t* pointset = bov_points_new((float[10][2]) {
	                                    {-1.0,  0.0},
	                                    {-0.8, -0.6},
	                                    {-0.7,  0.6},
	                                    {-0.5,  0.0},
	                                    {-0.2, -0.4},
	                                    { 0.0,  0.8},
	                                    { 0.3,  0.0},
	                                    { 0.5, -0.6},
	                                    { 0.7,  0.6},
	                                    { 0.0, -0.9}
	                                }, 10, GL_STATIC_DRAW);
	bov_points_set_color(pointset, (float[4]) {0.05, 0.1, 0.2, 0.6});

	bov_order_t* order = bov_order_new((GLuint[10]) {4, 3, 6, 9, 1, 0, 2, 5, 8, 7},
	                           10, GL_STATIC_DRAW);

	unsigned long frameCount = 0;
	while(!bov_window_should_close(window)) {
		double wtime = bov_window_get_time(window);

		bov_text_draw(window, common_text);

		switch( (unsigned) wtime / 4 % 8) {
		case 0:
			bov_text_draw(window, bov_points_draw_label);
			bov_points_draw(window, pointset, 0, BOV_TILL_END);
			break;
		case 1:
			bov_text_draw(window, bov_lines_draw_label);
			bov_lines_draw(window, pointset, 0, BOV_TILL_END);
			break;
		case 2:
			bov_text_draw(window, bov_line_strip_draw_label);
			bov_line_strip_draw(window, pointset, 0, BOV_TILL_END);
			break;
		case 3:
			bov_text_draw(window, bov_line_loop_draw_label);
			bov_line_loop_draw(window, pointset, 0, BOV_TILL_END);
			break;
		case 4:
			bov_text_draw(window, bov_curve_draw_label);
			bov_curve_draw(window, pointset, 0, BOV_TILL_END);
			break;
		case 5:
			bov_text_draw(window, bov_lines_draw_with_order_label);
			bov_lines_draw_with_order(window, pointset, order, 0, BOV_TILL_END);
			break;
		case 6:
			bov_text_draw(window, bov_line_strip_draw_with_order_label);
			bov_line_strip_draw_with_order(window, pointset, order, 0, BOV_TILL_END);
			break;
		case 7:
			bov_text_draw(window, bov_line_loop_draw_with_order_label);
			bov_line_loop_draw_with_order(window, pointset, order, 0, BOV_TILL_END);
			break;
		}

		bov_window_update(window);

		frameCount++;
	}

	printf("Ended correctly - %.2f second, %lu frames, %.2f fps\n",
	       bov_window_get_time(window),
	       frameCount,
	       frameCount / bov_window_get_time(window));

	bov_text_delete(bov_points_draw_label);
	bov_text_delete(bov_lines_draw_label);
	bov_text_delete(bov_line_strip_draw_label);
	bov_text_delete(bov_line_loop_draw_label);
	bov_text_delete(bov_curve_draw_label);
	bov_text_delete(bov_lines_draw_with_order_label);
	bov_text_delete(bov_line_strip_draw_with_order_label);
	bov_text_delete(bov_line_loop_draw_with_order_label);

	bov_points_delete(pointset);

	bov_order_delete(order);

	bov_window_delete(window);

	return EXIT_SUCCESS;
}
