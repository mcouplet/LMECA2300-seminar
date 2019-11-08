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

layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

flat out vec2 center;

void main() {
    vec2 center_scaling = localScale*(space_type==2 ? 2.0/resolution : scale); // 2/resolution: the size of one screen pixel
    vec2 width_scaling = localScale*(space_type!=0 ? 2.0/resolution : scale);
    vec2 center_translate = space_type==2 ? localPos : scale*(localPos + translate);

    center = gl_in[0].gl_Position.xy*center_scaling + center_translate;
    vec2 upRight = center + width*width_scaling;
    vec2 downLeft = center - width*width_scaling;
    vec2 upLeft = vec2(downLeft.x, upRight.y);
    vec2 downRight = vec2(upRight.x, downLeft.y);

    gl_Position = vec4(upLeft, 0.0, 1.0);
    EmitVertex();
    gl_Position = vec4(downLeft, 0.0, 1.0);
    EmitVertex();
    gl_Position = vec4(upRight, 0.0, 1.0);
    EmitVertex();
    gl_Position = vec4(downRight, 0.0, 1.0);
    EmitVertex();
    EndPrimitive();
}