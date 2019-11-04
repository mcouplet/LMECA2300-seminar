#include "draw_tools.h"
/***************************************
 *   FONT DEFINITION                   *
 ***************************************/
#include "font.h"

/***************************************
 *   SHADERS                           *
 ***************************************/
static GLchar points_vert[] = (GLchar[]) {
#include "points_vert.h"
};

static GLchar points_geom[] = (GLchar[]) {
#include "points_geom.h"
};

static GLchar points_frag[] = (GLchar[]) {
#include "points_frag.h"
};

static GLchar lines_geom[] = (GLchar[]) {
#include "lines_geom.h"
};

static GLchar lines_frag[] = (GLchar[]) {
#include "lines_frag.h"
};

static GLchar curve_geom[] = (GLchar[]) {
#include "curve_geom.h"
};

static GLchar text_vert[] = (GLchar[]) {
#include "text_vert.h"
};

static GLchar text_frag[] = (GLchar[]) {
#include "text_frag.h"
};


#define POS_LOCATION 0
#define TEX_LOCATION 1

#define TEXT_PROGRAM_INDEX     0
#define POINTS_PROGRAM_INDEX   1
#define LINES_PROGRAM_INDEX    2
#define CURVE_PROGRAM_INDEX    3
#define TRIANGLE_PROGRAM_INDEX 4
#define QUAD_PROGRAM_INDEX     5



/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Struct def
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
struct order_struct{
    GLuint ebo;
    size_t eboCapacity;
    size_t eboLen;
};

/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Errors
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
void error_log(int errorCode, const char *fmt, ...)
{
    if(errorCode>=OUT_OF_MEM_ERROR)
        perror("Current system error message: ");

    va_list vl;
    fprintf(stderr, "=X= Error %d: ", errorCode);
    va_start(vl,fmt);
    vfprintf(stderr, fmt, vl);
    va_end(vl);
    fputc('\n', stderr);
}

/* used to check the result of a memory allocation, exit on failure */
#define CHECK_MALLOC(ptr) if((ptr)==NULL) { \
        ERROR_LOG(OUT_OF_MEM_ERROR, "Memory allocation failed"); exit(EXIT_FAILURE); }


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Loading shaders
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
static int checkivError(GLuint object, GLenum pname,
                        PFNGLGETPROGRAMIVPROC glGetiv,
                        PFNGLGETPROGRAMINFOLOGPROC glGetInfoLog)
{
    GLint compile_status = GL_TRUE;
    glGetiv(object, pname, &compile_status);
    if(compile_status != GL_TRUE)
    {
        GLint logsize;
        glGetiv(object, GL_INFO_LOG_LENGTH, &logsize);

        char *log = malloc(logsize + 1);
        CHECK_MALLOC(log);

        glGetInfoLog(object, logsize, &logsize, log);
        if(pname==GL_LINK_STATUS)
            ERROR_LOG(SHADER_ERROR, "%s\t-> link operation failed", log);
        else if(pname==GL_COMPILE_STATUS)
            ERROR_LOG(SHADER_ERROR, "%s\t-> compile operation failed",log);
        else
            ERROR_LOG(SHADER_ERROR, "%s\t-> Unknown object operation failed",log);
        free(log);
        return SHADER_ERROR;
    }
    return 0;
}


static GLuint LoadShader(GLsizei count, const GLchar* shaderSource[], const GLint *lengths, const char* shaderName, GLenum shaderType)
{    
    GLuint shader = glCreateShader(shaderType);
    if(shader==0) {
        ERROR_LOG(SHADER_ERROR, "shader '%s' creation failed", shaderName);
        return 0;
    }

    glShaderSource(shader, count, shaderSource, lengths);
    glCompileShader(shader);

    if(checkivError(shader,GL_COMPILE_STATUS,glGetShaderiv,glGetShaderInfoLog)){
        const char* shaderType_desc = "unknown-type";
        if (shaderType==GL_VERTEX_SHADER)
            shaderType_desc = "vertex";
        else if (shaderType==GL_FRAGMENT_SHADER)
            shaderType_desc = "fragment";
        else if (shaderType==GL_GEOMETRY_SHADER)
            shaderType_desc = "geometry";
        ERROR_LOG(SHADER_ERROR, "%s shader: %s compilation failed", shaderType_desc, shaderName);
        return 0;
    }

    return shader;
}



static int program_init(window_t* window, int program_index, int n, ...)
{
    GLuint shaderProgram =  window->program[program_index];
    int i;
    va_list vl;

    va_start(vl,n);
    for (i=0; i<n; i++)
    {
        glAttachShader(shaderProgram, va_arg(vl,GLuint));
    }
    va_end(vl);

    glBindFragDataLocation(shaderProgram, 0, "outColor");
    glLinkProgram(shaderProgram);

    if(checkivError(shaderProgram,GL_LINK_STATUS,glGetProgramiv,glGetProgramInfoLog)){
        ERROR_LOG(SHADER_ERROR, "shader program %d creation failed", program_index);
        return 1;
    }
    
    return 0;
}

static GLuint create_texture(GLsizei width, GLsizei height, const GLvoid * data, GLenum format, GLint filterParam, GLint wrapParam)
{
    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);

    glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filterParam);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filterParam);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapParam);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapParam);

    // glBindTexture(GL_TEXTURE_2D, 0);

    return texture;
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Rasterizers
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
static void text_rasterizer_init(window_t* window) {
    window->texture = create_texture(font.tex_width, font.tex_height, font.tex_data, GL_RED, GL_LINEAR, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, window->texture);

    GLuint program = window->program[TEXT_PROGRAM_INDEX];
    glUseProgram(program);

    unsigned worldBlockIndex = glGetUniformBlockIndex(program, "worldBlock");
    unsigned objectBlockIndex = glGetUniformBlockIndex(program, "objectBlock");
    glUniformBlockBinding(program, worldBlockIndex, 0);
    glUniformBlockBinding(program, objectBlockIndex, 1);

    window->texture_sloc = glGetUniformLocation(program, "fontTex");
}

static void points_rasterizer_init(GLuint program) {
    glUseProgram(program);

    unsigned worldBlockIndex = glGetUniformBlockIndex(program, "worldBlock");
    unsigned objectBlockIndex = glGetUniformBlockIndex(program, "objectBlock");
    glUniformBlockBinding(program, worldBlockIndex, 0);
    glUniformBlockBinding(program, objectBlockIndex, 1);
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %   mouse_button_callback
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
static void mouse_button_callback(GLFWwindow* self, int button, int action, int mod)
{
    window_t* window = (window_t*) glfwGetWindowUserPointer(self);

    if(button==GLFW_MOUSE_BUTTON_LEFT){
        if(action==GLFW_PRESS){
            glfwSetCursor(self, window->leftClickCursor);
            window->clickTime[0] = window->wtime;
        }
        else{
            glfwSetCursor(self, NULL);
            window->clickTime[0] = -window->wtime;
        }
    }
    else if(button==GLFW_MOUSE_BUTTON_RIGHT){
        if(action==GLFW_PRESS){
            glfwSetInputMode(self, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
            window->clickTime[1] = window->wtime;
        }
        else{
            glfwSetInputMode(self, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
            window->clickTime[1] = -window->wtime;
        }
    }
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %   mouse_button_callback
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
static void scroll_callback(GLFWwindow* self, double x, double y) {
    window_t* window = (window_t*) glfwGetWindowUserPointer(self);

    if(y>0.0) {
        window->param.scale[0] *= 1 + .1*y;
        window->param.scale[1] *= 1 + .1*y;
    }
    else if(y<0.0) {
        window->param.scale[0] /= 1-.1*y;
        window->param.scale[1] /= 1-.1*y;
    }
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %   cursor_pos_callback
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
static void cursor_pos_callback(GLFWwindow* self, double x, double y)
{
    window_t* window = (window_t*) glfwGetWindowUserPointer(self);

    // TODO: modify this to use the scale
    float newX = (2.0*x/window->size[0]-1.0)/window->param.scale[0];
    float newY = (2.0*(1.0 - y/window->size[1])-1.0)/window->param.scale[1];

    if(window->clickTime[0]>0.0 || window->clickTime[1]>0.0) {
        window->param.translate[0] += (newX - window->cursorPos[0]);
        window->param.translate[1] += (newY - window->cursorPos[1]);
    }

    window->cursorPos[0] = newX;
    window->cursorPos[1] = newY;
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %      error_callback
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
static void error_callback(int errorCode, const char* description)
{
    ERROR_LOG(errorCode, "%s", description);
    exit(EXIT_FAILURE);
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %       key_callback
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
static void key_callback(GLFWwindow* self, int key, int scancode, int action, int mods)
{
    static unsigned screenshot_nbr = 0;
    char screenshot_name[64] = "screenshot";
    window_t* window = (window_t*) glfwGetWindowUserPointer(self);
    if(action==GLFW_PRESS || action==GLFW_REPEAT){
        switch(key){
            case GLFW_KEY_ESCAPE:
               glfwSetWindowShouldClose(self,GL_TRUE);
               break;
            case GLFW_KEY_SPACE:
                if(window->running==0) {
                    window->running = 1;
                }
                else {
                    window->running = 0;
                }
                break;
            case GLFW_KEY_P:
                snprintf(screenshot_name + 10, 54, "%u.ppm", screenshot_nbr++);
                window_screenshot(window, screenshot_name);
                break;
        }
    }
    if (key==GLFW_KEY_ESCAPE)
       glfwSetWindowShouldClose(self,GL_TRUE);
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  FrameBuffer_callback
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
static void framebuffer_size_callback(GLFWwindow* self, int width, int height)
{
    window_t* window = (window_t*) glfwGetWindowUserPointer(self);
    window->param.res[0] = width;
    window->param.res[1] = height;
    glViewport(0, 0, width, height);

    if(width>height) {
        window->param.scale[0] = window->param.scale[1]*height/width;
    }
    else {
        window->param.scale[1] = window->param.scale[0]*width/height;
    }
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Window_callback
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
static void window_size_callback(GLFWwindow* self, int width, int height)
{
    window_t* window = (window_t*) glfwGetWindowUserPointer(self);
    window->size[0] = width;
    window->size[1] = height;
}




/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Window object
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
static void window_OpenGL_init(window_t* window) {
    // Create draw and object UBO
    glGenBuffers(2, window->ubo);
    glBindBuffer(GL_UNIFORM_BUFFER, window->ubo[0]);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(world_param_t), NULL, GL_STATIC_DRAW);
    glBindBufferRange(GL_UNIFORM_BUFFER, 0, window->ubo[0], 0, sizeof(world_param_t));

    glBindBuffer(GL_UNIFORM_BUFFER, window->ubo[1]);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(object_param_t), NULL, GL_DYNAMIC_DRAW);
    glBindBufferRange(GL_UNIFORM_BUFFER, 1, window->ubo[1], 0, sizeof(object_param_t));
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    // compile shaders
    {
        GLuint fontVS, fontFS;
        if((fontVS = LoadShader(1,
                               (const GLchar* []) {text_vert},
                               (GLint[]){ sizeof(text_vert)-1 },
                               "font_vert.glsl",
                               GL_VERTEX_SHADER))==0)
            goto shader_error;
        if((fontFS = LoadShader(1,
                               (const GLchar* []) {text_frag},
                               (GLint[]){ sizeof(text_frag)-1 },
                               "font_frag.glsl",
                               GL_FRAGMENT_SHADER))==0)
            goto shader_error;

        window->program[TEXT_PROGRAM_INDEX] = glCreateProgram();

        // Specify the layout of the vertex data
        glBindAttribLocation(window->program[TEXT_PROGRAM_INDEX], POS_LOCATION, "pos");
        glBindAttribLocation(window->program[TEXT_PROGRAM_INDEX], TEX_LOCATION, "tex");

        if(program_init(window, TEXT_PROGRAM_INDEX, 2, fontVS, fontFS))
            goto shader_error;
        glDetachShader(window->program[TEXT_PROGRAM_INDEX],fontFS);
        glDetachShader(window->program[TEXT_PROGRAM_INDEX],fontVS);
    }

    {
        GLuint pointsVS, pointsGS, linesGS, curveGS, pointsFS, linesFS;

        if((pointsVS = LoadShader(1,
                                  (const GLchar* []) {points_vert},
                                  (const GLint[]) {sizeof(points_vert)-1},
                                  "points_vert.glsl",
                                  GL_VERTEX_SHADER))==0)
            goto shader_error;
        if((pointsGS = LoadShader(1,
                                  (const GLchar* []) {points_geom},
                                  (const GLint[]) {sizeof(points_geom)-1},
                                  "points_geo.glsl",
                                  GL_GEOMETRY_SHADER))==0)
            goto shader_error;
        if((linesGS = LoadShader(1,
                                  (const GLchar* []) {lines_geom},
                                  (const GLint[]) {sizeof(lines_geom)-1},
                                  "lines_geo.glsl",
                                  GL_GEOMETRY_SHADER))==0)
            goto shader_error;
        if((curveGS = LoadShader(1,
                                  (const GLchar* []) {curve_geom},
                                  (const GLint[]) {sizeof(curve_geom)-1},
                                  "curve_geo.glsl",
                                  GL_GEOMETRY_SHADER))==0)
            goto shader_error;
        if((pointsFS = LoadShader(1,
                                  (const GLchar* []) {points_frag},
                                  (const GLint[]) {sizeof(points_frag)-1},
                                  "points_frag.glsl",
                                  GL_FRAGMENT_SHADER))==0)
            goto shader_error;
        if((linesFS = LoadShader(1,
                                  (const GLchar* []) {lines_frag},
                                  (const GLint[]) {sizeof(lines_frag)-1},
                                  "lines_frag.glsl",
                                  GL_FRAGMENT_SHADER))==0)
            goto shader_error;

        window->program[POINTS_PROGRAM_INDEX] = glCreateProgram();
        window->program[LINES_PROGRAM_INDEX] = glCreateProgram();
        window->program[CURVE_PROGRAM_INDEX] = glCreateProgram();

        // Specify the layout of the vertex data
        glBindAttribLocation(window->program[POINTS_PROGRAM_INDEX], POS_LOCATION, "pos");
        glBindAttribLocation(window->program[LINES_PROGRAM_INDEX], POS_LOCATION, "pos");
        glBindAttribLocation(window->program[CURVE_PROGRAM_INDEX], POS_LOCATION, "pos");

        if(program_init(window, POINTS_PROGRAM_INDEX,
                        3, pointsVS, pointsGS, pointsFS))
            goto shader_error;
        if(program_init(window, LINES_PROGRAM_INDEX,
                        3, pointsVS, linesGS, linesFS))
            goto shader_error;
        if(program_init(window, CURVE_PROGRAM_INDEX,
                        3, pointsVS, curveGS, linesFS))
            goto shader_error;

        glDetachShader(window->program[POINTS_PROGRAM_INDEX],pointsVS);
        glDetachShader(window->program[POINTS_PROGRAM_INDEX],pointsGS);
        glDetachShader(window->program[POINTS_PROGRAM_INDEX],pointsFS);
        glDetachShader(window->program[LINES_PROGRAM_INDEX],pointsVS);
        glDetachShader(window->program[LINES_PROGRAM_INDEX],linesGS);
        glDetachShader(window->program[LINES_PROGRAM_INDEX],linesFS);
        glDetachShader(window->program[CURVE_PROGRAM_INDEX],pointsVS);
        glDetachShader(window->program[CURVE_PROGRAM_INDEX],curveGS);
        glDetachShader(window->program[CURVE_PROGRAM_INDEX],linesFS);
    }

    // Create all rasterizers
    text_rasterizer_init(window);
    points_rasterizer_init(window->program[POINTS_PROGRAM_INDEX]);
    points_rasterizer_init(window->program[LINES_PROGRAM_INDEX]);
    points_rasterizer_init(window->program[CURVE_PROGRAM_INDEX]);

    glUseProgram(window->program[LINES_PROGRAM_INDEX]);
    window->last_program = LINES_PROGRAM_INDEX;

    return;

shader_error:
    ERROR_LOG(SHADER_ERROR, "check your driver and OpenGL capabilities");
    exit(EXIT_FAILURE);
}


window_t* window_new(int width, int height, const char* win_name)
{
    window_t* window = malloc(sizeof(window_t));
    CHECK_MALLOC(window);

    glfwSetErrorCallback(error_callback);

    if (!glfwInit())
        exit(EXIT_FAILURE);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_SAMPLES, 4); // multisampling (MSAA) x4
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    if(width==0 || height==0){
        const GLFWvidmode *mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
        window->size[0] = mode->width;
        window->size[1] = mode->height;
        window->self = glfwCreateWindow(mode->width, mode->height, win_name, glfwGetPrimaryMonitor(), NULL); // fullscreen
    }
    else{
        if(width<0) {
            glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE);
            width = -width;
        }
        window->size[0] = width;

        if(height<0) {
            glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
            height = -height;
        }
        window->size[1] = height;
        
        window->self = glfwCreateWindow(width, height, win_name, NULL, NULL);
    }

    if (!window->self) {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window->self);
    glfwSetWindowUserPointer(window->self,window);

    // load opengl functions with GLAD
    if(!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress))
        ERROR_LOG(GLAD_ERROR, "Failed to initialize OpenGL context");

    glfwGetFramebufferSize(window->self, &width, &height);

    window->param.scale[0] = window->param.scale[1] = 1.0;
    framebuffer_size_callback(window->self, width, height);
    window->param.translate[0] = window->param.translate[1] = 0.0;
    // window->param.rotation = 0.0;
    window->wtime = 0.0;


    // glfwSetWindowCloseCallback(window->self,close_callback);
    glfwSetKeyCallback(window->self, key_callback);
    glfwSetFramebufferSizeCallback(window->self, framebuffer_size_callback);
    glfwSetWindowSizeCallback(window->self, window_size_callback);
    glfwSetMouseButtonCallback(window->self, mouse_button_callback);
    glfwSetCursorPosCallback(window->self, cursor_pos_callback);
    glfwSetScrollCallback(window->self, scroll_callback);

    //glfwSetInputMode(window->self,GLFW_CURSOR,GLFW_CURSOR_HIDDEN);



    glfwSwapInterval(1);            // vsync
    // clear screen to white
    glClearColor(1.f,1.f,1.f,1.00);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_CULL_FACE);     // cull face that are not counterclockwise automatically
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    glfwSetTime(0.0);              // set the time to 0, avoid overflow when converting to GLfloat...

    {
        double cursorX = 0.0, cursorY = 0.0;
        glfwGetCursorPos(window->self, &cursorX, &cursorY);
        window->cursorPos[0] = 2.0*cursorX/window->size[0]-1.0;
        window->cursorPos[1] = 2.0*(1.0 - cursorY/window->size[1])-1.0;
        window->clickTime[0] = 0.0;
        window->clickTime[1] = 0.0;
    }

    window->running = 1;
    window->leftClickCursor = glfwCreateStandardCursor(GLFW_HRESIZE_CURSOR);

    window_OpenGL_init(window);

    return window;
}


void window_update(window_t* window) {
    // Swap front and back buffers
    glfwSwapBuffers(window->self);

    if(window->running) {
        // Poll for and process events
        glfwPollEvents();
        window->wtime = glfwGetTime();
    }
    else {
        glfwWaitEvents();
        glfwSetTime(window->wtime);
    }

    // update the world ubo (we update every frame, even if it isn't changed...)
    glBindBuffer(GL_UNIFORM_BUFFER, window->ubo[0]);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(world_param_t), &window->param);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    glClear(GL_COLOR_BUFFER_BIT);
}

void window_update_and_wait_events(window_t* window) {
    int state = window->running;
    window->running = 0;
    window_update(window);
    window->running = state;
}


void window_delete(window_t* window)
{
    glDeleteTextures(1, &window->texture);
    glDeleteProgram(window->program[TEXT_PROGRAM_INDEX]);
    glDeleteProgram(window->program[POINTS_PROGRAM_INDEX]);
    glDeleteProgram(window->program[LINES_PROGRAM_INDEX]);
    glDeleteProgram(window->program[CURVE_PROGRAM_INDEX]);
    glDeleteBuffers(2, window->ubo);
    glfwDestroyWindow(window->self);
    glfwTerminate();
    free(window);
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Order object
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
order_t* order_new(GLuint* elements, size_t n, GLenum usage) {
    order_t* order = malloc(sizeof(order_t));
    CHECK_MALLOC(order);

    order->eboLen = elements==NULL ? 0 : n;

    glGenBuffers(1, &order->ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, order->ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, n*sizeof(GLuint), elements, GL_STATIC_DRAW);
    order->eboCapacity = n;
    return order;
}


order_t* order_update(order_t* order, GLuint* elements, size_t n) {
    order->eboLen = elements==NULL ? 0 : n;

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, order->ebo);
    if(order->eboLen > order->eboCapacity){
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, n*sizeof(GLuint), elements, GL_DYNAMIC_DRAW);
        order->eboCapacity = order->eboLen;
    }
    else if(elements!=NULL){
        glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, n*sizeof(GLuint), elements);
    }

    return order;
}

order_t* order_partial_update(order_t* order, GLuint* elements, size_t start, size_t end, size_t newN) {
    if(elements==NULL) {
        ERROR_LOG(PARAMETER_ERROR, "Cannot do a partial update whith a NULL pointer for array of elements");
        return NULL;
    }

    if(end > newN)
        newN = end;

    if(newN > order->eboCapacity) {
        ERROR_LOG(PARAMETER_ERROR, "Cannot do a partial update when the new size is bigger than the capacity of the buffer");
        return NULL;
    }

    order->eboLen = newN;

    if(start>end)
        return order;

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, order->ebo);
    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, start, (end-start)*sizeof(GLuint), elements);

    return order;
}


void order_delete(order_t* order) {
    glDeleteBuffers(1, &order->ebo);
    free(order);
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Text Object
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
static size_t fill_text_data(float* data, const unsigned char* string, size_t len)
{
    float pen_x = 0;
    float pen_y = 0;
    float x,y,w,h;
    // float direction = 1.0;

    size_t num = 0;
    for (size_t i=0; i<len; i++) {
        texture_glyph_t *glyph = font.glyphs + string[i];

        switch(glyph->codepoint){
            // up, semi-up, semi-down...  x01 - 06, e1 - e6, c0 -df, f0 - ff  could mean something
            case '\f':  // back to the beginning
                pen_y = 0;
            case '\n':  // newline
            
            case '\r':  // carriage return
                pen_x = 0;
            case ' ':   // space
            case '\a':  // goes up :-)
            case '\v':  // vertical space
            case '\b':  // back the size of a space
                break;
            case '\t':  // tab
                pen_x = ((int) (pen_x/font.size + 2)/2)*2 * font.size;
                break;
            default:
                x = (pen_x + glyph->offset_x)/font.size;
                y = (pen_y + glyph->offset_y)/font.size;

                w = glyph->width/font.size;
                h = glyph->height/font.size;

                data[24*num+0] = x;
                data[24*num+1] = y;
                data[24*num+2] = glyph->s0;
                data[24*num+3] = glyph->t0;
                data[24*num+4] = x;
                data[24*num+5] = y-h;
                data[24*num+6] = glyph->s0;
                data[24*num+7] = glyph->t1;
                data[24*num+8] = x+w;
                data[24*num+9] = y-h;
                data[24*num+10] = glyph->s1;
                data[24*num+11] = glyph->t1;
                data[24*num+12] = x;
                data[24*num+13] = y;
                data[24*num+14] = glyph->s0;
                data[24*num+15] = glyph->t0;
                data[24*num+16] = x+w;
                data[24*num+17] = y-h;
                data[24*num+18] = glyph->s1;
                data[24*num+19] = glyph->t1;
                data[24*num+20] = x + w;
                data[24*num+21] = y;
                data[24*num+22] = glyph->s1;
                data[24*num+23] = glyph->t0;

                num++;
                break;

        }
        pen_x += glyph->advance_x;
        pen_y += glyph->advance_y;
    }

    return num;
}

// size (in characters) to alloc, and usage (either GL_STATIC_DRAW or GL_DYNAMIC_DRAW)
text_t* text_new(unsigned char* string, GLenum usage){
    text_t* text = malloc(sizeof(text_t));
    CHECK_MALLOC(text);

    text->param = (object_param_t){
                    0.0, 0.0, 0.0, 1.0, // color
                    1.0 ,1.0, 1.0, 1.0, // outlineColor
                    0.0 ,0.0,           // other
                    0.0, 0.0,           // localPos
                    0.05, 0.05,         // localScale
                    0.0,                // width
                    -1.0,               // outlineWidth
                    // 0.0,                // rotation
                    NORMAL_SPACE};

    // Create Vertex Array Object
    glGenVertexArrays(1, &text->vao);
    glBindVertexArray(text->vao);

    // Vertex Buffer Object
    glGenBuffers(1, &text->vbo);
    glBindBuffer(GL_ARRAY_BUFFER, text->vbo);

    // specify the layout of the data
    glVertexAttribPointer(POS_LOCATION, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), 0);
    glVertexAttribPointer(TEX_LOCATION, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), (void*)(2*sizeof(float)));
    glEnableVertexAttribArray(POS_LOCATION);
    glEnableVertexAttribArray(TEX_LOCATION);

    if(string!=NULL) {
        text->dataCapacity = strlen(string);

        // data for text contain... 4 vertex per letter
        // for each letter, we must have
        // 2 screen coordinates      4f
        // 2 texture coordinates     4f

        text->data = malloc(24*text->dataCapacity*sizeof(float));
        CHECK_MALLOC(text->data);

        // text->string = string;

        text->vboLen = fill_text_data(text->data, string, text->dataCapacity);
        glBufferData(GL_ARRAY_BUFFER, 24*text->vboLen*sizeof(float), text->data, usage);
        text->vboCapacity = text->vboLen;
    }
    else {
        text->vboLen = 0;
        text->vboCapacity = 0;
    }

    // glBindBuffer(GL_ARRAY_BUFFER, 0);
    // glBindVertexArray(0);

    return text;
}

void text_delete(text_t* text){
    glDeleteBuffers(1, &text->vbo);
    glDeleteVertexArrays(1, &text->vao);
    free(text->data);
    free(text);
}

text_t* text_update(text_t* text, unsigned char* string){
    // see if the length is not longer than the original string
    size_t newLen = strlen(string);

    if(newLen > text->dataCapacity){
        free(text->data);
        text->data = malloc(24*newLen*sizeof(float));
        CHECK_MALLOC(text->data);
        text->dataCapacity = newLen;
    }

    
    // text->string = string;
    text->vboLen = fill_text_data(text->data, string, newLen);

    glBindBuffer(GL_ARRAY_BUFFER, text->vbo);
    if(text->vboLen > text->vboCapacity){
        glBufferData(GL_ARRAY_BUFFER, 24*text->vboLen*sizeof(float), text->data, GL_DYNAMIC_DRAW);
        text->vboCapacity = text->vboLen;
    }
    else{
        glBufferSubData(GL_ARRAY_BUFFER, 0, 24*text->vboLen*sizeof(float), text->data);
    }

    return text;
}


void text_draw(window_t* window, text_t* text){
    if(text->vboLen==0)
        return;

    if(window->last_program!=TEXT_PROGRAM_INDEX) {
        glUseProgram(window->program[TEXT_PROGRAM_INDEX]);
        window->last_program = TEXT_PROGRAM_INDEX;
    }

    // update the object ubo
    glBindBuffer(GL_UNIFORM_BUFFER, window->ubo[1]);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(object_param_t), &text->param);
    // glBindBuffer(GL_UNIFORM_BUFFER, 0);

    glUniform1i(window->texture_sloc, 0);
    glBindVertexArray(text->vao);
    glDrawArrays(GL_TRIANGLES, 0, text->vboLen*6);
    // glBindVertexArray(0);
}



/*%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Points Object
 %%%%%%%%%%%%%%%%%%%%%%%%%*/
points_t* points_new(float* coords, size_t n, GLenum usage) {
    points_t* points = malloc(sizeof(points_t));
    CHECK_MALLOC(points);

    points->vboLen = coords==NULL ? 0: n;

    points->param = (object_param_t){
                    0.0, 0.0, 0.0, 1.0, // color
                    1.0 ,1.0, 1.0, 1.0, // outlineColor
                    0.0 ,0.0,           // other
                    0.0, 0.0,           // localPos
                    1.0, 1.0,           // localScale
                    0.025,              // width
                    -1.0,               // outlineWidth
                    // 0.0,                // rotation
                    NORMAL_SPACE};

    // Create Vertex Array Object
    glGenVertexArrays(1, &points->vao);
    glBindVertexArray(points->vao);

    // Vertex Buffer Object
    glGenBuffers(1, &points->vbo);
    glBindBuffer(GL_ARRAY_BUFFER, points->vbo);

    // specify the layout of the data
    glVertexAttribPointer(POS_LOCATION, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(POS_LOCATION);

    glBufferData(GL_ARRAY_BUFFER, 2*n*sizeof(float), coords, usage);
    points->vboCapacity = n;

    // glBindBuffer(GL_ARRAY_BUFFER, 0);
    // glBindVertexArray(0);

    return points;
}

points_t* points_update(points_t* points, float* coords, size_t n) {
    points->vboLen = coords==NULL ? 0: n;

    glBindBuffer(GL_ARRAY_BUFFER, points->vbo);
    if(n > points->vboCapacity){
        glBufferData(GL_ARRAY_BUFFER, 2*n*sizeof(float), coords, GL_DYNAMIC_DRAW);
        points->vboCapacity = n;
    }
    else if(coords!=NULL){
        glBufferSubData(GL_ARRAY_BUFFER, 0, 2*n*sizeof(float), coords);
    }
    // glBindBuffer(GL_ARRAY_BUFFER, 0);

    return points;
}

points_t* points_partial_update(points_t* points, float* coords, size_t start, size_t end, size_t newN) {
    if(coords==NULL) {
        ERROR_LOG(PARAMETER_ERROR, "Cannot do a partial update whith a NULL pointer for array of coordinates");
        return NULL;
    }

    if(end > newN)
        newN = end;

    if(newN > points->vboCapacity) {
        ERROR_LOG(PARAMETER_ERROR, "Cannot do a partial update when the new size is bigger than the capacity of the buffer");
        return NULL;
    }

    points->vboLen = newN;

    if(start>end)
        return points;

    glBindBuffer(GL_ARRAY_BUFFER, points->vbo);
    glBufferSubData(GL_ARRAY_BUFFER, start, (end-start)*2*sizeof(float), coords);
    // glBindBuffer(GL_ARRAY_BUFFER, 0);

    return points;
}

void points_delete(points_t* points){
    glDeleteBuffers(1, &points->vbo);
    glDeleteVertexArrays(1, &points->vao);
    free(points);
}


static inline void switch_rasterizer_with_mode(window_t* window, GLenum mode)
{
    int index;

    switch(mode) {
        case GL_POINTS:
            index = POINTS_PROGRAM_INDEX;
            break;
        case GL_LINES :
        case GL_LINE_STRIP :
        case GL_LINE_LOOP :
            index = LINES_PROGRAM_INDEX;
            break;
        case GL_LINE_STRIP_ADJACENCY :
            index = CURVE_PROGRAM_INDEX;
            break;
        case GL_LINES_ADJACENCY:
            index = QUAD_PROGRAM_INDEX;
            break;
        case GL_TRIANGLES:
        case GL_TRIANGLE_STRIP:
            index = TRIANGLE_PROGRAM_INDEX;
            break;
        default:
            ERROR_LOG(PARAMETER_ERROR, "this function is private but was given erroneous arguments anyway");
            exit(EXIT_FAILURE);
    }

    if(window->last_program!=index) {
        glUseProgram(window->program[index]);
        window->last_program = index;
    }
}

void points_draw_aux(window_t* window, points_t* points, GLenum mode) {
    if(points->vboLen==0)
        return;

    switch_rasterizer_with_mode(window, mode);

    // update the object ubo
    glBindBuffer(GL_UNIFORM_BUFFER, window->ubo[1]);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(object_param_t), &points->param);
    // glBindBuffer(GL_UNIFORM_BUFFER, 0);

    glBindVertexArray(points->vao);
    glDrawArrays(mode, 0, points->vboLen);
    // glBindVertexArray(0);
}

void points_draw_with_order_aux(window_t* window, points_t* points, order_t* order, GLenum mode) {
    if(points->vboLen==0 || order->eboLen==0)
        return;

    switch_rasterizer_with_mode(window, mode);

    // update the object ubo
    glBindBuffer(GL_UNIFORM_BUFFER, window->ubo[1]);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(object_param_t), &points->param);
    // glBindBuffer(GL_UNIFORM_BUFFER, 0);

    glBindVertexArray(points->vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, order->ebo);
    glDrawElements(mode, order->eboLen, GL_UNSIGNED_INT, 0);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    // glBindVertexArray(0);
}

void window_screenshot(window_t* window, char* filename) {
    static unsigned char* data = NULL;

    if(window==NULL || filename==NULL) {
        free(data);
        return;
    }

    int width = window->param.res[0];
    int height = window->param.res[1];

    int rowsize = (width*3+3)/4*4; // width of a row in byte
    data = realloc(data, rowsize*height); // width is rounded to multiple of four

    FILE *pFile;

    pFile=fopen(filename, "wb");
    if(pFile==NULL)
        exit(1);

    glReadBuffer( GL_FRONT );
    glReadPixels( 0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, data );

    // Write header
    fprintf(pFile, "P6\n%d %d\n255\n", width, height);

    // Write pixel data
    for(int i=height-1; i>=0; i--){
        unsigned char *row = data + i*rowsize;
        fwrite(row, 1, width*3, pFile);
    }

    // Close file
    fclose(pFile);
}