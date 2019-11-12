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

// gets a position in world space,
// return the same position in sceen space 
vec2 getScreenPos(vec2 pos) {
    vec2 resRatio = min(resolution.x, resolution.y)/resolution;
    if(space_type==0) {
        // classical case
        return resRatio*zoom*(localScale*pos + localPos + translate);
    }
    else if(space_type==1) {
        // no scaling of the height is applied
        return resRatio*(zoom*(localPos+translate)+localScale*pos);
    }
    else {
        // everything is given in pixel, from the bottom left corner :-)
        return 2.0*(localScale*pos + localPos)/resolution - 1.0;
    }

    // getScreenScaledVec(worldPos) + getScreenTranslation();
}

void main() {
    vec2 a = gl_in[0].gl_Position.xy;
    vec2 b = gl_in[1].gl_Position.xy;

    vec2 ba = b - a;
    lba = length(ba);
    vec2 v = (ba)/lba; // direction of this segment
    vec2 n = vec2(-v.y, v.x);

    float minRes = min(resolution.x, resolution.y);
    vec2 resRatio = minRes/resolution;

    // screenPos = scaling*worldPos + translation

    // localScale should not affect the width and outlineWidth
    // therefore widthScaling!=scaling

    // pixelSize is simply 2.0/resolution / widthScaling
    vec2 scaling;
    vec2 widthScaling;
    vec2 translation;
    if(space_type==0) {
        // classical case
        widthScaling = resRatio*zoom;
        scaling = resRatio*zoom*localScale;
        translation = resRatio*zoom*(localPos + translate);
        pixelSize = 2.0/(minRes*zoom);
    }
    else if(space_type==1) {
        widthScaling = resRatio; // same as 0 but no zoom
        scaling = resRatio*localScale; // same as 0 but no zoom
        translation = resRatio*zoom*(localPos + translate); // same as 0
        pixelSize = 2.0/minRes;
    }
    else {
        widthScaling = 2.0/resolution;
        scaling = localScale*widthScaling;
        translation = localPos*widthScaling - 1.0;
        pixelSize = 1.0;
    }

    pRect = vec2(-width);
    gl_Position = vec4(scaling*a - widthScaling*width*(v+n) + translation, 0.0, 1.0);
    EmitVertex();
    pRect = vec2(lba + width, -width);
    gl_Position = vec4(scaling*b + widthScaling*width*(v-n) + translation, 0.0, 1.0);
    EmitVertex();
    pRect = vec2(-width, width);
    gl_Position = vec4(scaling*a + widthScaling*width*(n-v) + translation, 0.0, 1.0);
    EmitVertex();
    pRect = vec2(lba + width, width);
    gl_Position = vec4(scaling*b + widthScaling*width*(v+n) + translation, 0.0, 1.0);
    EmitVertex();
    EndPrimitive();
}