#version 150 core

layout(std140) uniform objectBlock
{
    vec4 fillColor;
    vec4 outlineColor;
    vec2 otherParam; // x: pointiness y: ?
    vec2 localPos;
    vec2 localScale;
    float width;
    float outlineWidth;
    // float rotation;
    int space_type; // 0: normal sizes, 1: size in pixels, 2: size in pixels without translation
};

layout(std140) uniform worldBlock
{
    vec2 resolution;
    vec2 scale;
    vec2 translate;
    // float rotation;
    // float wtime;
};

flat in vec2 p1_screen;
flat in vec2 p2_screen;

out vec4 outColor;

#define antialiasing 4.0 // number of pixel for the antialiasing

// return shortest vector from point p to segment ab, given ba and pa
vec2 dist2Segment(vec2 ba, vec2 pa) {
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return pa - ba*h;
}

void main( void ) {
    vec2 pixel = 2.0*gl_FragCoord.xy/resolution - 1.0;

    vec2 center_scaling = localScale*(space_type==2 ? 2.0/resolution : scale);// 2/resolution: the size of one screen pixel
    vec2 width_scaling = localScale*(space_type!=0 ? 2.0/resolution : scale);

    vec2 d = dist2Segment((p2_screen-p1_screen)/center_scaling, (pixel - p1_screen)/width_scaling);
    float len = length(d);

    // compute the length of the smoothstep in pixel
    float smoothLength = len==0.0 ? 0.0 : antialiasing*len/length(width_scaling*d*resolution);

    // float a=abs(line_dist);
    float opacity = smoothstep(width,width-smoothLength,len);

    float mu = smoothstep(width-outlineWidth, width-outlineWidth-smoothLength, len);
    outColor = mix(outlineColor, fillColor, mu); // at 0: completely outlineColor, at1: completely fillColor
    outColor.a *= opacity;
}