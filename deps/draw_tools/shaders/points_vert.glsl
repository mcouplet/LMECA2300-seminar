#version 150 core

layout (std140) uniform objectBlock
{
    vec4 fillColor;
    vec4 outlineColor;
    vec2 localPos;
    vec2 localScale;
    float width;
    float marker;
    float outlineWidth;
    // float rotation;
    int space_type; // 0: normal sizes, 1: unzoomable, 2: unmodifable pixel size
};

layout (std140) uniform worldBlock
{
    vec2 resolution;
    vec2 scale;
    vec2 translate;
    // float rotation;
};

in vec2 pos;

void main()
{
    if(space_type<2) {
        // classical case
        gl_Position = vec4(scale*localScale*(localPos + translate), 0.0, 1.0);
    }
    else if(space_type==1) {
        // no scaling of the height is applied
        gl_Position = vec4(scale*localScale*(localPos + translate), 0.0, 1.0);
    }
    else {
        // everything is given in pixel, from the bottom left corner :-)
        gl_Position = vec4(localScale*localPos/resolution - 1.0, 0.0, 1.0);
    }
}