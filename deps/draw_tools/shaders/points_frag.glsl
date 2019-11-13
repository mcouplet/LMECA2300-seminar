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

in vec2 pSquare;
flat in float pixelSize;

out vec4 outColor;

#define INV_SQRT_2 0.70710678118

void main( void ) {
    vec2 pixel = 2.0*gl_FragCoord.xy/resolution - 1.0;

    vec2 width_scaling = localScale*(space_type!=0 ? 2.0/resolution : scale);

    vec2 circle = (center - pixel)/width_scaling;
    float amplitude_circle = length(circle);
    vec2 direction_circle = normalize(circle);

    float xS = sign(circle.x);
    float yS = sign(circle.y);
    float xL = xS*circle.x; // or abs(circle.x)
    float yL = yS*circle.y; // or abs(circle.y)

    float amplitude_square = max(xL,yL);
    vec2 direction_square = xL>yL ? vec2(xS, 0.0) : vec2(0.0, yS);

    float amplitude_rhombus = (INV_SQRT_2*(xL+yL)+(1.-INV_SQRT_2)*width);
    vec2 direction_rhombus = vec2(xS, yS)*INV_SQRT_2;

    float amplitude_emptyCross = min(xL,yL)+width;
    vec2 direction_emptyCross = xL<yL ? vec2(xS, 0.0) : vec2(0.0, yS);

    float amplitude_invCircle = 2.0*width - length(circle - width*vec2(xS, yS));
    vec2 direction_invCircle = normalize(circle - width*vec2(xS, yS));

    // rhombus - (circle-rhombus)
    // float amplitude_weird1 = 2.0*amplitude_rhombus - amplitude_circle;
    // vec2 direction_weird1 = normalize(2.0*direction_rhombus - direction_circle);

    // // circle - (square - circle)
    // float amplitude_weird2 = 2.0*amplitude_circle - amplitude_square;
    // vec2 direction_weird2 = normalize(2.0*direction_circle - direction_square);

    vec2 direction;
    float amplitude;

    /* possibilities:
       rhombus - square
       rhombus - circle
       rhombus - emptyCross
       rhombus - invCircle
       square-circle
       
          0        1           2           3          4        5       6=>0
       square - rhombus - emptyCross - invCircle - rhombus - circle - square
    */

    float markerMod6 = mod(marker, 6);

    if(markerMod6 < 1.0) {
        direction = normalize(mix(direction_square, direction_rhombus, markerMod6 - 0.0));
        amplitude = mix(amplitude_square, amplitude_rhombus, markerMod6 - 0.0);
    }
    else if(markerMod6 < 2.0) {
        direction = normalize(mix(direction_rhombus, direction_emptyCross, markerMod6 - 1.0));
        amplitude = mix(amplitude_rhombus, amplitude_emptyCross, markerMod6 - 1.0);
    }
    else if(markerMod6 < 3.0) {
        direction = normalize(mix(direction_emptyCross, direction_invCircle, markerMod6 - 2.0));
        amplitude = mix(amplitude_emptyCross, amplitude_invCircle, markerMod6 - 2.0);
    }
    else if(markerMod6 < 4.0) {
        direction = normalize(mix(direction_invCircle, direction_rhombus, markerMod6 - 3.0));
        amplitude = mix(amplitude_invCircle, amplitude_rhombus, markerMod6 - 3.0);
    }
    else if(markerMod6 < 5.0) {
        direction = normalize(mix(direction_rhombus, direction_circle, markerMod6 - 4.0));
        amplitude = mix(amplitude_rhombus, amplitude_circle, markerMod6 - 4.0);
    }
    else {
        direction = normalize(mix(direction_circle, direction_square, markerMod6 - 5.0));
        amplitude = mix(amplitude_circle, amplitude_square, markerMod6 - 5.0);
    }

    // compute the length of the smoothstep in pixel
    float smoothLength = amplitude==0.0 ? 0.0 : antialiasing/length(width_scaling*direction*resolution);

    // float a=abs(line_dist);
    float opacity = smoothstep(width, width-smoothLength,amplitude);

    float mu = smoothstep(width-outlineWidth, width-outlineWidth-smoothLength, amplitude);
    outColor = mix(outlineColor, fillColor, mu); // at 0: completely outlineColor, at1: completely fillColor
    // outColor = mix(vec4(0.0,0.0,0.0,1.0), outColor, opacity);
    outColor.a *= opacity;

    // outColor = vec4(amplitude/width);
}