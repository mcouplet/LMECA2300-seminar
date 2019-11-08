#ifndef __DRAW_TOOLS_H__
#define __DRAW_TOOLS_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

/* opaque structures */
typedef struct window_struct window_t;
typedef struct order_struct order_t;
typedef struct text_struct text_t;
typedef struct points_struct points_t;


/***
 *      _____         _      _  _                   _____  _____ 
 *     |  __ \       | |    | |(_)           /\    |  __ \|_   _|
 *     | |__) |_   _ | |__  | | _   ___     /  \   | |__) | | |  
 *     |  ___/| | | || '_ \ | || | / __|   / /\ \  |  ___/  | |  
 *     | |    | |_| || |_) || || || (__   / ____ \ | |     _| |_ 
 *     |_|     \__,_||_.__/ |_||_| \___| /_/    \_\|_|    |_____|
 *                                                               
 *                                                               
 */

/* a space_type_t can be given to the functions text_set_space_type() and
 * points_set_space_type() */
typedef enum {
    USUAL_SPACE = 0,      // you can zoom and translate the space in which
                          // the object is embeded

    UNZOOMABLE_SPACE = 1, // zooming will not change the size of the object

    PIXEL_SPACE = 2,      // unzoomable, unmovable and positon and scaling must
                          // be given in pixel coordinates
} space_type_t;


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Window
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
/* creating a window is probably the first thing a program must do.
 * width is either:
 *     - the width of the window in pixel
 *     - 0 for full screen
 *     - a negative number for a maximized window
 * height is either:
 *     - the height of the window in pixel
 *     - 0 for full screen
 *     - negative height for a fixed size window
 * win_name is the title of the window
 * return a window object
 */
window_t* window_new(int width, int height, const char* win_name);

/* Once something was drawn to the framebuffer, the window must be updated to
 * display the content of the framebuffer on the screen. Updating the window
 * also handles the mouse and keyboard inputs, thus you must redraw your scene
 * frequently even if nothing changed in order to get smooth input handling.
 */
void window_update(window_t* window);

/* same as window_update but wait for input events before clearing the screen */
void window_update_and_wait_events(window_t* window);

/* delete the window properly */
void window_delete(window_t* window);


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Text
 %%%%%%%%%%%%%%%%%%%%%%%%%*/

/* 
 * Create a text object from a string. string a null-terminated C string usage
 * is either:
 *    - GL_STATIC_DRAW if you don't intend to change the content of the text
 *      object (the string of characters)
 *    - GL_DYNAMIC_DRAW if you intend to change it with the text_update function
 */
text_t* text_new(unsigned char* string, GLenum usage);

/* change the content of the text object */
text_t* text_update(text_t* text, unsigned char* string);

/* draw a text object, using a text rasterizer (note: a text rasterizer can be
 * used to draw multiple text object)*/
void text_draw(window_t* window, text_t* text);

/* delete text properly */
void text_delete(text_t* text);


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Points
 %%%%%%%%%%%%%%%%%%%%%%%%%*/

/* Create a points object
 * coord is either:
 *    - an array of coordinates, interleaved : x1, y1, x2, y2, ... xn, yn
 *         => n is the number of points and the maximum capacity
 *    - NULL, in which case a points object is created with a maximum capacity
 *      of n points usage is either:
 *    - GL_STATIC_DRAW if you don't intend to update the coordinates of the
 *      points regularly
 *    - GL_DYNAMIC_DRAW if you intend to change coordinates regularly with the
 *      points_update or points_partial_update function
 */
points_t* points_new(GLfloat* coords, GLsizei n, GLenum usage);

/* change the content of the points object (the maximum capacity can be increased)*/
points_t* points_update(points_t* points, GLfloat* coords, GLsizei n);

/* change the content of the points object, but only from start to end (the
 * maximum capacity cannot be increased)*/
points_t* points_partial_update(points_t* points, GLfloat* coords, GLsizei start, GLsizei end, GLsizei newN);

/* draw points markers to the window */
static inline void points_draw(window_t* window, points_t* pts);
static inline void points_draw_with_order(window_t* window, points_t* pts, order_t* order);

/* draw lines between pairs of points, (p1->p2) (p3->p4) (p5->p6)... */
static inline void lines_draw(window_t* window, points_t* pts);
static inline void lines_draw_with_order(window_t* window, points_t* pts, order_t* order);

/* draw lines that connect each points p1->p2->p3->p4->p5->p6...->pn */
static inline void line_strip_draw(window_t* window, points_t* pts);
static inline void line_strip_draw_with_order(window_t* window, points_t* pts, order_t* order);

/* draw lines that connect each points and end with the first
 * p1->p2->p3->p4->p5->p6...->pn->p1*/
static inline void line_loop_draw(window_t* window, points_t* pts);
static inline void line_loop_draw_with_order(window_t* window, points_t* pts, order_t* order);

/* draw a single line that connect each points p1->p2->p3->p4->p5->p6...->pn
 * The difference with line_strip_draw can really be seen with an outline or
 * transparency */
static inline void curve_draw(window_t* window, points_t* pts);
static inline void curve_draw_with_order(window_t* window, points_t* pts, order_t* order) ;

/* delete a point object */
void points_delete(points_t* points);


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Order
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
/* 
 * Create an order object, that enables to draw only certain points, or to draw
 * lines that connect different points elements is either:
 *     - an array of indices (the indices of the point that we will draw) => n
 *       is the number of indices
 *     - NULL, in which case an order object is created with a maximum capacity
 *       of n indices usage is either:
 *    - GL_STATIC_DRAW if you don't intend to update the content of the order
 *      object regularly
 *    - GL_DYNAMIC_DRAW if you intend to change it regularly with the
 *      order_update or order_partial_update function
 */
order_t* order_new(GLuint* elements, GLsizei n, GLenum usage);

/* change the content of the order object (the maximum capacity can be
 * increased) */
order_t* order_update(order_t* order, GLuint* elements, GLsizei n, GLenum usage);

/* change the content of the points object, but only from start to end (the
 * maximum capacity cannot be increased) */
order_t* order_partial_update(order_t* order, GLuint* elements, GLsizei start, GLsizei end, GLsizei newN);

/* delete an order object */
void order_delete(order_t* order);


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Window parameters
 %%%%%%%%%%%%%%%%%%%%%%%%%*/

/* get the time in second associated with the window. The time is updated via
 * the window_update() and window_update_and_wait_events() functions. The time
 * is stopped if the user press the space bar. */
static inline double window_get_time(window_t* window);

/* get the resolution in pixel of the window */
static inline GLfloat window_get_xres(window_t* window);
static inline GLfloat window_get_yres(window_t* window);

/* tells if the window should close (because the user decided it) */
static inline int window_should_close(window_t* window);
static inline void window_set_color(window_t* window, GLfloat rgba[4]);

/* translates the content of the window (similar to dragging with the mouse */
static inline void window_translate(window_t* window, GLfloat posx, GLfloat posy);

/* scale the content of the window (similar to zooming) */
static inline void window_scale(window_t* window, GLfloat scale);

/* take a screenshot and save it as a PPM with the name 'filename' */
void window_screenshot(window_t* window, char* filename);

/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  text parameters
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
/* set the position of the text on the screen.
 * Note that coordinates must be given in pixels if the space type is set to
 * PIXEL_SPACE*/
static inline void text_set_pos(text_t* text, GLfloat pos_x, GLfloat pos_y);

/* set the height of the text object on the screen.
 * Note that this height must be given in pixels if
 * the space type is set to PIXEL_SPACE */
static inline void text_set_height(text_t* text, GLfloat height);

// not implemented yet
// static inline void text_set_rotation(text_t* text, GLfloat rotation);

/* sets the text color to a 4 channel color (red, gree, blue, alpha) */
static inline void text_set_color(text_t* text, GLfloat rgba[4]);

/* boldness is a parameter, usually between -1 and 1, that define
 * how massive the font should be. -1 means light, 1 means bold */ 
static inline void text_set_boldness(text_t* text, GLfloat width);

static inline void text_set_outline_color(text_t* text, GLfloat rgba[4]);

/* for text, the outline width must be between 0 and 1
 * 0 means no outline, 1 means the outline might meet in the middle of the character */
static inline void text_set_outline_width(text_t* text, GLfloat width);

/* the outline can be slightly shifted in some direction.
 * More specifically, if you have a normal vector to the font 'n',
 * and a shift vector 's', and an outline width 'w', then we will have
 * a new outline width 'w = w + dot(s,n)'.
 * For example, if we have an outline width of 0.2, and a shift vector of
 * (0.1, 0.0), the outline will be 0.3 on the right of each character and
 * 0.1 on the left of each character
 *
 * Note: font with an outline shift can become very bad when zoom in.
 * actually, it uses the gradient of the sdf based font to calculate
 * the normal, and this gradient approaches 0 when zoomed in...
 */
static inline void text_set_outline_shift(text_t* text, GLfloat shift_x, GLfloat shift_y);

/* for more detail on what this function does, see:
 * - the structure `space_type_t`
 * - the function `text_set_pos()` and `text_set_height()` */
static inline void text_set_space_type(text_t* text, space_type_t space_type);


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  points parameters
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
/* set the position of the points on the screen.
 * Note that coordinates must be given in pixels if the space type is set to
* PIXEL_SPACE*/
static inline void points_set_pos(points_t* points, GLfloat pos_x, GLfloat pos_y);

/* you can give a local scaling factor to your set of points in both x and y
 * directions.*/
static inline void points_scale(points_t* points, GLfloat scaleX, GLfloat scaleY);

// not implemented yet
// static inline void points_set_rotation(points_t* points, GLfloat rotation);

/* set the width of the marker/line/curve/... associated with the points */
static inline void points_set_width(points_t* points, GLfloat width);

/* sets color of the marker/line/curve/... associated with the points to a 4 
 * channel color (red, gree, blue, alpha) */
static inline void points_set_color(points_t* points, GLfloat rgba[4]);

// not implemented yet
// static inline void points_set_pointiness(points_t* points, GLfloat pointiness);


static inline void points_set_outline_color(points_t* points, GLfloat rgba[4]);
static inline void points_set_outline_width(points_t* points, GLfloat width);
static inline void points_set_marker(points_t* points, GLfloat marker);
static inline void points_set_space_type(points_t* points, space_type_t space_type);

static inline void program_wait_events();


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %      error_log
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
/* some error code in addition to those defined by GLFW (www.glfw.org/docs/latest/group__errors.html) */
#define PARAMETER_ERROR   (0x00020001)
#define SHADER_ERROR      (0x00020002)
#define GLAD_ERROR        (0x00020003)

#define OUT_OF_MEM_ERROR  (0x00030001)
#define IO_ERROR          (0x00030002)

/* display an error with ERROR_LOG, for example :
 * if(parameter_not_valid) {
 *     ERROR_LOG(PARAMETER_ERROR, "parameter XX must be positive !");
 *     return NULL;
 */
#ifndef NDEBUG
#define ERROR_LOG(errorCode, fmt, ...) do{ error_log(errorCode, fmt, ## __VA_ARGS__ ); \
        fprintf(stderr, "\t(in function %s, line %d)\n", __func__, __LINE__); }while(0)
#else
#define ERROR_LOG(errorCode,fmt, ...) error_log(errorCode, fmt, ## __VA_ARGS__ )
#endif


/***********************************************************************************/














/*
 *      _____   _____   _____ __      __    _______  ______ 
 *     |  __ \ |  __ \ |_   _|\ \    / //\ |__   __||  ____|
 *     | |__) || |__) |  | |   \ \  / //  \   | |   | |__   
 *     |  ___/ |  _  /   | |    \ \/ // /\ \  | |   |  __|  
 *     | |     | | \ \  _| |_    \  // ____ \ | |   | |____ 
 *     |_|     |_|  \_\|_____|    \//_/    \_\|_|   |______|
 *
 *
 *       __| | ___  _ __ |/| |_  | | ___   ___ | | __ | |__   ___| | _____      __
 *      / _` |/ _ \| '_ \  | __| | |/ _ \ / _ \| |/ / | '_ \ / _ \ |/ _ \ \ /\ / /
 *     | (_| | (_) | | | | | |_  | | (_) | (_) |   <  | |_) |  __/ | (_) \ V  V /
 *      \__,_|\___/|_| |_|  \__| |_|\___/ \___/|_|\_\ |_.__/ \___|_|\___/ \_/\_/
 */

typedef struct world_param_struct world_param_t;
typedef struct object_param_struct object_param_t;

struct world_param_struct{
    GLfloat res[2];
    GLfloat scale[2];
    GLfloat translate[2];
    // GLfloat rotation;
    // GLfloat wtime;
};

struct object_param_struct {
    GLfloat fillColor[4];
    GLfloat outlineColor[4];
    GLfloat other[2];
    GLfloat localPos[2];
    GLfloat localScale[2];
    GLfloat width;
    GLfloat outlineWidth;
    // GLfloat rotation;
    GLint space_type;
};


struct window_struct
{
    GLFWwindow* self;
    GLFWcursor* leftClickCursor;  // cursor for the left click

    world_param_t param;

    int size[2];         // size of window (might not be in pixel => !=framebuffer)
    double cursorPos[2]; // new position of the cursor in screen coordinates
    double clickTime[2]; // time of the click (left and right) or -1 if released
    double wtime;

    GLuint ubo[2];       // one ubo for a world_param_t, one ubo for an object_param_t

    int last_program;
    GLuint program[6];
    GLuint texture;
    int texture_sloc;

    int running;         // running or paused
};

struct text_struct {
    GLuint vao;
    GLuint vbo;
    GLsizei vboCapacity; // capacity of the vbo
    GLsizei vboLen;      // number of letter in vbo
    // const unsigned char* string;
    GLsizei dataCapacity; // length of the string
    GLfloat* data;
    object_param_t param;
};

struct points_struct{
    GLuint vao;
    GLuint vbo;
    GLsizei vboCapacity; // numer of points in vbo
    GLsizei vboLen;
    object_param_t param;
};


static inline void program_wait_events() {
    glfwWaitEvents();
}

static inline double window_get_time(window_t* window){
    return window->wtime;
}

static inline GLfloat window_get_xres(window_t* window){
    return window->param.res[0];
}

static inline GLfloat window_get_yres(window_t* window){
    return window->param.res[1];
}


static inline int window_should_close(window_t* window){
    return glfwWindowShouldClose(window->self);
}

static inline void window_set_color(window_t* window, GLfloat rgba[4]) {
    glClearColor(rgba[0], rgba[1], rgba[2], rgba[3]);
}

static inline void window_translate(window_t* window, GLfloat posx, GLfloat posy) {
    window->param.translate[0] += posx;
    window->param.translate[1] += posy;
}

static inline void window_scale(window_t* window, GLfloat scale) {
    window->param.scale[0] *= scale;
    window->param.scale[1] *= scale;
}


static inline void text_set_pos(text_t* text, GLfloat pos_x, GLfloat pos_y) {
    text->param.localPos[0] = pos_x;
    text->param.localPos[1] = pos_y;
}

static inline void text_set_height(text_t* text, GLfloat height) {
    text->param.localScale[0] = height;
    text->param.localScale[1] = height;
}

// static inline void text_set_rotation(text_t* text, GLfloat rotation) {
//     text->param.rotation = rotation;
// }

static inline void text_set_color(text_t* text, GLfloat rgba[4]) {
    for (int i=0; i<4; i++)
        text->param.fillColor[i] = rgba[i];
}

static inline void text_set_boldness(text_t* text, GLfloat boldness) {
    text->param.width = boldness;
}

static inline void text_set_outline_color(text_t* text, GLfloat rgba[4]) {
    for (int i=0; i<4; i++)
        text->param.outlineColor[i] = rgba[i];
}

static inline void text_set_outline_width(text_t* text, GLfloat width) {
    text->param.outlineWidth = width;
}

static inline void text_set_outline_shift(text_t* text, GLfloat shift_x, GLfloat shift_y) {
    text->param.other[0] = shift_x;
    text->param.other[1] = shift_y;
}

static inline void text_set_space_type(text_t* text, space_type_t space_type) {
    text->param.space_type = space_type;
}



static inline void points_set_pos(points_t* points, GLfloat pos_x, GLfloat pos_y) {
    points->param.localPos[0] = pos_x;
    points->param.localPos[1] = pos_y;
}

static inline void points_scale(points_t* points, GLfloat scaleX, GLfloat scaleY) {
    points->param.localScale[0] = scaleX;
    points->param.localScale[1] = scaleY;
}

// static inline void points_set_rotation(points_t* points, GLfloat rotation) {
//     points->param.rotation = rotation;
// }

static inline void points_set_width(points_t* points, GLfloat width) {
    points->param.width = width;
}

static inline void points_set_color(points_t* points, GLfloat rgba[4]) {
    for (int i=0; i<4; i++)
        points->param.fillColor[i] = rgba[i];
}

// a number between 0 and 1
// static inline void points_set_pointiness(points_t* points, GLfloat pointiness) {
//     points->param.other[0] = pointiness;
// }

static inline void points_set_outline_color(points_t* points, GLfloat rgba[4]) {
    for (int i=0; i<4; i++)
        points->param.outlineColor[i] = rgba[i];
}

static inline void points_set_outline_width(points_t* points, GLfloat width) {
    points->param.outlineWidth = width;
}

static inline void points_set_marker(points_t* points, GLfloat marker) {
    points->param.other[1] = marker;
}

static inline void points_set_space_type(points_t* points, space_type_t space_type) {
    points->param.space_type = space_type;
}

void points_draw_aux(window_t* window, points_t* pts, GLenum mode);
void points_draw_with_order_aux(window_t* window, points_t* pts, order_t* order, GLenum mode);

static inline void points_draw(window_t* window, points_t* pts) {
    points_draw_aux(window, pts, GL_POINTS);
}

static inline void points_draw_with_order(window_t* window, points_t* pts, order_t* order) {
    points_draw_with_order_aux(window, pts, order, GL_POINTS);
}

static inline void curve_draw(window_t* window, points_t* pts) {
    points_draw_aux(window, pts, GL_LINE_STRIP_ADJACENCY);
}

static inline void curve_draw_with_order(window_t* window, points_t* pts, order_t* order) {
    points_draw_with_order_aux(window, pts, order, GL_LINE_STRIP_ADJACENCY);
}

static inline void lines_draw(window_t* window, points_t* pts) {
    points_draw_aux(window, pts, GL_LINES);
}

static inline void lines_draw_with_order(window_t* window, points_t* pts, order_t* order) {
    points_draw_with_order_aux(window, pts, order, GL_LINES);
}

static inline void line_strip_draw(window_t* window, points_t* pts) {
    points_draw_aux(window, pts, GL_LINE_STRIP);
}

static inline void line_strip_draw_with_order(window_t* window, points_t* pts, order_t* order) {
    points_draw_with_order_aux(window, pts, order, GL_LINE_STRIP);
}

static inline void line_loop_draw(window_t* window, points_t* pts) {
    points_draw_aux(window, pts, GL_LINE_LOOP);
}

static inline void line_loop_draw_with_order(window_t* window, points_t* pts, order_t* order) {
    points_draw_with_order_aux(window, pts, order, GL_LINE_LOOP);
}

static inline void triangles_draw(window_t* window, points_t* pts) {
    points_draw_aux(window, pts, GL_TRIANGLES);
}

static inline void triangles_draw_with_order(window_t* window, points_t* pts, order_t* order) {
    points_draw_with_order_aux(window, pts, order, GL_TRIANGLES);
}

static inline void triangle_strip_draw(window_t* window, points_t* pts) {
    points_draw_aux(window, pts, GL_TRIANGLE_STRIP);
}

static inline void triangle_strip_draw_with_order(window_t* window, points_t* pts, order_t* order) {
    points_draw_with_order_aux(window, pts, order, GL_TRIANGLE_STRIP);
}

static inline void quad_draw(window_t* window, points_t* pts) {
    points_draw_aux(window, pts, GL_LINES_ADJACENCY);
}

static inline void quad_draw_with_order(window_t* window, points_t* pts, order_t* order) {
    points_draw_with_order_aux(window, pts, order, GL_LINES_ADJACENCY);
}

#ifdef __GNUC__
void error_log(int errorCode, const char *fmt, ...)__attribute__((format (printf, 2, 3)));
#else
void error_log(int errorCode, const char *fmt, ...);
#endif

#endif

