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
    float vl;
    // we should do the same as for the points :-)
    // if(shape>1.0) {
        vl = length(v); // circle end
    // }
    // else if(shape<1.0){
        // vl = max(abs(v.x), abs(v.y)); // square end
    // }
    // else {
        // vl = abs(v.x) + abs(v.y); // triangle end
    // }

    vec2 w = vec2(width-outlineWidth, width);
    vec2 alpha = smoothstep(w+pixelSize, w-pixelSize, vec2(vl));
    outColor = mix(outlineColor, fillColor, alpha.x); // at 0: completely outlineColor, at1: completely fillColor
    outColor.a *= alpha.y;
}