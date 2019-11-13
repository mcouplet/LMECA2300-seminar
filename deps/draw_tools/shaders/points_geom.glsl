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

layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

out vec2 pSquare;
flat out float pixelSize;

void main() {
    vec2 p = gl_in[0].gl_Position.xy*localScale;

    float minRes = min(resolution.x, resolution.y);
    vec2 resRatio = minRes/resolution;

    // screenPos = scaling*worldPos + translation

    // localScale should not affect the width and outlineWidth
    // therefore widthScaling!=scaling

    // pixelSize is simply 2.0/resolution / scaling
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

    vec2 center = p*scaling + translation;
    vec2 upRight = center + width*scaling;
    vec2 downLeft = center - width*scaling;
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