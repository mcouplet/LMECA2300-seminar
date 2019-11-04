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
    // float wtime;
};

layout(lines) in;
layout(triangle_strip, max_vertices = 4) out;

flat out vec2 p1_screen;
flat out vec2 p2_screen;

vec2 perpendicular(vec2 v) {
    return vec2(-v.y, v.x);
}

void main() {
    vec2 p1 = gl_in[0].gl_Position.xy;
    vec2 p2 = gl_in[1].gl_Position.xy;

    vec2 center_scaling = localScale*(space_type==2 ? 2.0/resolution : scale);// 2/resolution: the size of one screen pixel
    vec2 width_scaling = localScale*(space_type!=0 ? 2.0/resolution : scale);
    vec2 center_translate = space_type==2 ? localPos : scale*(localPos + translate);

    vec2 v1=normalize(p2 - p1); // direction of this segment
    vec2 n1 = perpendicular(v1);

    p1_screen = p1*center_scaling+center_translate;
    p2_screen = p2*center_scaling+center_translate;
    gl_Position = vec4(p1_screen - width_scaling*width*(v1-n1), 0.0, 1.0);
    EmitVertex();
    gl_Position = vec4(p1_screen - width_scaling*width*(v1+n1), 0.0, 1.0);
    EmitVertex();
    gl_Position = vec4(p2_screen + width_scaling*width*(v1+n1), 0.0, 1.0);
    EmitVertex();
    gl_Position = vec4(p2_screen + width_scaling*width*(v1-n1), 0.0, 1.0);
    EmitVertex();
    EndPrimitive();
}