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

void main( void ) {
    float shape = mod(marker, 3.0f);
    vec2 v = vec2(pRect.x - clamp(pRect.x, 0.0, lba), pRect.y);
    vec2 sdf = length(v) - vec2(width, width-outlineWidth+step(outlineWidth, 0.0)); // circles

    vec2 alpha = smoothstep(pixelSize, -pixelSize, sdf);
    outColor = mix(outlineColor, fillColor, alpha.y); // at 0: completely outlineColor, at1: completely fillColor
    outColor.a *= alpha.x;
}