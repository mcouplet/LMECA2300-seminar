#version 150 core

layout (std140) uniform objectBlock
{
    vec4 fillColor;
    vec4 outlineColor;
    vec2 otherParam;
    vec2 localPos;
    vec2 localScale;
    float width;
    float outlineWidth;
    // float rotation;
    int space_type; // 0: normal sizes, 1: size in pixels, 2: size in pixels without translation
};

layout (std140) uniform worldBlock
{
    vec2 resolution;
    vec2 scale;
    vec2 translate;
    // float rotation;
};

out vec2 texCoord;

in vec2 pos; // this position is in pixel, starting from the top...
in vec2 tex;

void main()
{
    texCoord = tex;

    vec2 width_scaling = localScale*(space_type!=0 ? 2.0/resolution : scale);
    vec2 center_translate = space_type==2 ? localPos : scale*(localPos + translate);

    gl_Position = vec4(pos*width_scaling + center_translate, 0.0, 1.0); // scaling is already in translation
}