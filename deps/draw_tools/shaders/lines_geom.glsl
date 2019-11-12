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
    vec2 translate;
    float zoom;
    // float rotation;
};

layout(lines) in;
layout(triangle_strip, max_vertices = 4) out;

out vec2 pRect;
flat out float pixelSize;
flat out float lba;

void main() {
    vec2 a = gl_in[0].gl_Position.xy*localScale;
    vec2 b = gl_in[1].gl_Position.xy*localScale;

    vec2 ba = b - a;
    lba = length(ba);
    vec2 v = ba/lba;         // direction of this segment
    vec2 n = vec2(-v.y, v.x);// perpendicular direction

    float minRes = min(resolution.x, resolution.y);
    vec2 resRatio = minRes/resolution;

    // screenPos = scaling*worldPos + translation

    // localScale should not affect the width and outlineWidth
    // therefore widthScaling!=scaling

    // pixelSize is simply 2.0/resolution / widthScaling
    vec2 scaling;
    vec2 translation;
    if(space_type==0) {
        // classical case
        scaling = resRatio*zoom;
        translation = resRatio*zoom*(localPos + translate);
        pixelSize = 2.0/(minRes*zoom);
    }
    else if(space_type==1) {
        scaling = resRatio;      // same as 0 but no zoom
        translation = resRatio*(localPos + zoom*translate); // same as 0
        pixelSize = 2.0/minRes;
    }
    else /*if(space_type==2)*/{
        scaling = 2.0/resolution;
        translation = localPos*scaling - 1.0;
        pixelSize = 1.0;
    }

    vec2 aScreen = scaling*a + translation;
    vec2 bScreen = scaling*b + translation;
    vec2 wScreen = scaling*width*v;
    vec2 hScreen = scaling*width*n;

    pRect = vec2(-width);
    gl_Position = vec4(aScreen - wScreen - hScreen, 0.0, 1.0);
    EmitVertex();
    pRect = vec2(lba + width, -width);
    gl_Position = vec4(bScreen + wScreen - hScreen, 0.0, 1.0);
    EmitVertex();
    pRect = vec2(-width, width);
    gl_Position = vec4(aScreen - wScreen + hScreen, 0.0, 1.0);
    EmitVertex();
    pRect = vec2(lba + width, width);
    gl_Position = vec4(bScreen + wScreen + hScreen, 0.0, 1.0);
    EmitVertex();
    EndPrimitive();
}