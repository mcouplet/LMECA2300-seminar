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

in vec2 pRect;
flat in float lba;
flat in float pixelSize;

out vec4 outColor;

#define antialiasing (4.0*pixelSize) // number of pixel for the antialiasing

void main( void ) {
    vec2 v = vec2(pRect.x - clamp(pRect.x, 0.0, lba), pRect.y);
    float vl = length(v);

    float opacity = smoothstep(width, width-antialiasing,vl);
    float mu = smoothstep(width-outlineWidth, width-outlineWidth-antialiasing, vl);
    outColor = mix(outlineColor, fillColor, mu); // at 0: completely outlineColor, at1: completely fillColor
    outColor.a *= opacity;
}