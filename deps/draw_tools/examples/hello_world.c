#include "draw_tools.h"
#include <math.h>

#define PI 3.14159265358979323846

void nice_colormap(float color[4], float a)
{
   color[0] = sin(0.33 * a) * 0.3 + 0.7;
   color[1] = sin(0.23 * a + 2.0) * 0.3 + 0.7;
   color[2] = sin(0.17 * a + 5.0) * 0.3 + 0.6;
   color[3] = 1.0; // solid
}

int main()
{
    window_t* window = window_new(800, 800, "Tutorial 1");
    window_set_color(window, (GLfloat[]){0.5, 0.5, 0.5, 1.0});

    // hw is prefix for hello world :p
    text_t* hw_obj = text_new((unsigned char[]){"Hello World !"},
                              GL_STATIC_DRAW);

    float hw_color[4] = {1.0, 1.0, 1.0, 1.0};
    text_set_fontsize(hw_obj, 0.2);

    // a character here is 0.1 wide (it is 0.05 with no scaling)
    // centering "Hello World !": 13/2=6.5 characters => 0.65
    text_set_pos(hw_obj, (GLfloat[2]) {-0.65, 0.0});
    text_set_boldness(hw_obj, 0.4); // bold
    text_set_outline_width(hw_obj, 1.0); // big outline
    text_set_outline_color(hw_obj, (float[4]) {0, 0, 0, 1}); // black


    while(!window_should_close(window)) {
        nice_colormap(hw_color, window_get_time(window));
        text_set_color(hw_obj, hw_color);
        text_draw(window, hw_obj);

        window_update(window);
    }

    text_delete(hw_obj);
    window_delete(window);

    return EXIT_SUCCESS;
}