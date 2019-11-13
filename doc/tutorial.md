Doing animations in C with `draw_tools`
=======================================

This page teaches you how to use the `draw_tools.h` functions in order
to make a beautiful animation for your upcoming project.

 * [FAQ](#faq)
 * [Tutorial](#tutorial)
	- [Creating a window](#creating-a-window)
	- [A window has its own time](#a-window-has-its-own-time)
	- [Hello World](#hello-world)
	- [Hello real world](#hello-real-world)

## FAQ

**Why C ?**

The first question that you might ask is:

> Why is the Numerical Geometry project in C and not in Python, Java, 
  C++, Rust or Julia ?

You will realize that C is very transparent about memory allocation.
Every non-constant chunk of memory that you use must be allocated on
the heap via a `malloc()` call. Therefore, you will realize precisely
how much memory is consumed by your algorithm, and you will naturally
avoid unnecessary allocations. In addition, in 2D, it is very easy to
allocate the right amount of memory right from the beginning. Take a
Delaunay Triangulation algorithm for example. If you have *n* points,
you are going to end up with *2\*(n-2)* triangles. You can allocate
all the space needed right from the beginning. That transparency and
simplicity give C a clear advantage over languages with fancy (but
costly) features. Finally, a well-realized C program will be as fast
as any C++ or Rust.

**What is draw_tools ?**

It is a small and simple OpenGL library created by your assistants to
help you draw lines, text, and simple markers (little shapes to 'mark'
the position of a point). All functions are defined and briefly
explained in `draw_tools.h`. The most useful ones are described in the$
tutorial below. We figured it would be easier than using something
like CSFML or BGFX, and lighter than requiring the SDL.



## Tutorial


### Creating a window

First, we will create a simple window (an object of type `window_t`)
via a call to the function
```C
    window_t* window_new(int width, int height, const char* win_name);
```
where
 * **width** is the width of the window in pixels or 0 for full-screen
 * **height** is the height of the window in pixels or 0 for
   full-screen
 * **win_name** is a string which will be the name of your window

we will wait for 2 seconds using the function `time()` and
`difftime()` given in
[`<time.h>`](http://www.cplusplus.com/reference/ctime/)

and close the window (and deleted all structure related to the window)
with
```C
    void window_delete(window_t* window);
```


Therefore, in the following code, we create a window of size 800x800
titled "Tutorial 0", we wait 2 seconds, close the window and return.
```C
#include "draw_tools.h"
#include <time.h>

int main()
{
    window_t* window = window_new(800, 800, "Tutorial 0");

    time_t tic = time(NULL);

    while(difftime(time(NULL) , tic)< 2)
        ;// do nothing

    window_delete(window);

    return 0;
}
```


### A window has its own time

Before drawing things into a window, we have to understand that the
window object is actually much more than a window.
Indeed, when you will draw primitive shapes with `points_draw()`,
`text_draw()`... you will in fact draw to a texture in the memory of 
your computer which is called a framebuffer. Once you've filled the
framebuffer with what you want to draw, you can show it in the window
using:
```C
void window_update(window_t* window)
```
That function actually does a lot of things:

 - it swaps the framebuffer with the screen buffer.
 - it waits for a screen refresh if you have VSYNC. A refresh is the
   moment when the screen buffer is actually shown on your screen. It
   depends on the refreshing rate (framerate) of your screen (usually
   60 Hz)
 - it processes events (mouse inputs and keyboard inputs) that
   happened since the last call to `window_update()`.
 - it updates the size of the framebuffer in pixel according to the
   size of the window.
 - **it clears the framebuffer**
 - it changes the render position, change the scaling (zoom), pause
   the window timer, tell if the window should close etc. according to
   the event received.
 - it updates its own timer which doesn't increase if the window is
   paused


Now, let's update our previous example using the window timer:
```C
#include "draw_tools.h"

int main()
{
    window_t* window = window_new(800, 800, "mY wInDoW");

    while(window_get_time(window) < 2) {
        window_update(window);
    }

    window_delete(window);

    return EXIT_SUCCESS;
}
```

Because we use the timer of the window, **we can pause the time by
pressing the space bar**. If you do so before the 2 seconds delay, the
window won't close because the time will get stuck somewhere in
between 0 and 2 seconds. Pressing the space bar again will restart the
timer and the window will close.

You might also notice that the background color is now white, which is
the default color for the background of the framebuffer. In the
previous example, it stayed black because we were never showing the
content of the framebuffer to the window.

You can change the background color of the window using:
```C
void window_set_color(window_t* window, float rgba[4]);
```

The color `rgba` is defined by 4 floating-point values `float` in the
range 0.0 - 1.0.
 1. red
 2. green
 3. blue
 4. alpha (the transparency)

Note that the transparency is totally useless with the window
background color, because transparency only has an impact when things
are drawn behind our object.

---

**waiting for user to close the window**

Showing a window during 2 second is not very useful. You will
generally want to display the window until the user clicks the close
button, press the escape key or ALT+F4. Those events are captured by
the window when using `window_update()`. You can simply ask the window
if they happened by calling
```C
int window_should_close(window_t* window)
```
in this manner

```C
while(!window_should_close(window)) {
	// draw stuff here

	// ...

    window_update(window);
}
```



### Hello World

Time to render our first primitive!

Let's use the function `text_new()`, which returns a new text
object:
```C
text_t* text_new(unsigned char* string, GLenum usage);
```
 - **string** is the message that you want to display
 - **usage** is either GL_DYNAMIC_DRAW if you intend to change the
   content (the text) of the object regularly or GL_STATIC_DRAW if the
   text won't change.

Additionally, there are multiple parameters associated to a text
object that you can change: its position on the screen, its color, its
width, its outline color, its outline width, its scaling in x and y...
see `draw_tools.h` for more info.


To render our text object into the framebuffer, we have to call
```C
void text_draw(window_t* window, text_t* text);
```
Reminder: we will have to do this repeatedly because `window_update()`
clears the framebuffer.


Finally, when we don't need our text object anymore, we will call
```C
void text_delete(text_t* text);
```
to destroy it.



Using all your knowledge, you should be able to display "Hello world"
with some pretty effects, by varying the text object parameters.
```C
#include "draw_tools.h"
#include <math.h>

#define PI 3.14159265358979323846

void nice_colormap(float color[4], float a)
{
   color[0] = sin(.33*a)*0.3+0.7;
   color[1] = sin(.23*a+2.)*0.3+0.7;
   color[2] = sin(.17*a+5.)*0.3+0.6;
   color[3] = 1; // solid
}

int main()
{
    window_t* window = window_new(800,800, "Tutorial 1");
    window_set_color(window, (GLfloat[]){0.5, 0.5, 0.5, 1.0});

    // hw is prefix for hello world :p
    text_t* hw_obj = text_new((unsigned char[]){"Hello World !"},
                              GL_STATIC_DRAW);

    float hw_color[4] = {1.0, 1.0, 1.0, 1.0};
    text_set_height(hw_obj, 0.2);

    // a character here is 0.1 wide (it is 0.05 with no scaling)
    // centering "Hello World !": 13/2=6.5 characters => 0.65
    text_set_pos(hw_obj, (GLfloat[2]){-0.65, 0.0});
    text_set_boldness(hw_obj, 0.4); // bold
    text_set_outline_width(hw_obj, 1.0); // big outline
    text_set_outline_color(hw_obj, (float[4]) {0, 0, 0, 1}); // black


    while(!window_should_close(window)){
        nice_colormap(hw_color, window_get_time(window));
        text_set_color(hw_obj, hw_color);
        text_draw(window, hw_obj);

        window_update(window);
    }

    text_delete(hw_obj);
    window_delete(window);

    return EXIT_SUCCESS;
}
```


### Hello real world

It's pretty easy to display simple graphics, but displaying the
progress of an algorithm is far more difficult. Indeed, you are not
only programming an algorithm for your Numerical Geometry project,
**you are a video game developer now !**. Indeed, the window must be
responsive to user inputs, users must be able to move the scene with
their mouses, zooming in and out instantaneously. To do so, you should
never do something that last more than a hundredth of a second before
updating the window. In addition, we want your project to be pretty,
with nice animations. The speed of your animations should not depend
on the framerate of your screen. These are constraints that are almost
only seen in video game development.

---

**the next and last part of the tutorial should arrive soon**
