 /*************************************************************************
  * Draw_tools 0.1
  * A wrapper around OpenGL and GLFW (www.glfw.org) to draw simple 2D
  * graphics.
  *------------------------------------------------------------------------
  * Copyright (c) 2019-2020 Célestin Marot <marotcelestin@gmail.com>
  *
  * This software is provided 'as-is', without any express or implied
  * warranty. In no event will the authors be held liable for any damages
  * arising from the use of this software.
  *
  * Permission is granted to anyone to use this software for any purpose,
  * including commercial applications, and to alter it and redistribute it
  * freely, subject to the following restrictions:
  *
  * 1. The origin of this software must not be misrepresented; you must not
  *    claim that you wrote the original software. If you use this software
  *    in a product, an acknowledgment in the product documentation would
  *    be appreciated but is not required.
  *
  * 2. Altered source versions must be plainly marked as such, and must not
  *    be misrepresented as being the original software.
  *
  * 3. This notice may not be removed or altered from any source
  *    distribution.
  *
  *************************************************************************/

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

// see https://www.iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm

float ndot(vec2 a, vec2 b ) { return a.x*b.x - a.y*b.y; }
float sdRhombus( in vec2 p, in vec2 b ) 
{
    vec2 q = abs(p);
    float h = clamp((-2.0*ndot(q,b)+ndot(b,b))/dot(b,b),-1.0,1.0);
    float d = length( q - 0.5*b*vec2(1.0-h,1.0+h) );
    return d * sign( q.x*b.y + q.y*b.x - b.x*b.y );
}

// p is the coordinate, s is the size
float sdCross(vec2 p, vec2 s)
{
    p = abs(p); p = (p.y>p.x) ? p.yx : p.xy;
    vec2  q = p - s;
    float k = max(q.y,q.x);
    vec2  w = (k>0.0) ? q : vec2(s.y-p.x,-k);
    return sign(k)*length(max(w,0.0));
}

float sdBox(vec2 p, vec2 b)
{
    vec2 d = abs(p)-b;
    return length(max(d,vec2(0))) + min(max(d.x,d.y),0.0);
}

// signed distance to a n-star polygon with external angle en
float sdStar(vec2 p, float r, int n, float m) // m=[2,n]
{
    // these 4 lines can be precomputed for a given shape
    float an = 3.141593/float(n);
    float en = 3.141593/m;
    vec2  acs = vec2(cos(an),sin(an));
    vec2  ecs = vec2(cos(en),sin(en)); // ecs=vec2(0,1) and simplify, for regular polygon,

    // reduce to first sector
    float bn = mod(atan(p.x,p.y),2.0*an) - an;
    p = length(p)*vec2(cos(bn),abs(sin(bn)));

    // line sdf
    p -= r*acs;
    p += ecs*clamp( -dot(p,ecs), 0.0, r*acs.y/ecs.y);
    return length(p)*sign(p.x);
}

void main( void ) {
    float m = mod(marker, 25.);      // marker%23
    float f = fract(m);
    float fw = f*width;          // fictive width
    float delta = fw - width;    // the difference in distance

    int shape = int(m);

    float sdf;
    if(shape==0)
        sdf = length(pSquare) - width;
    else if(shape==1)
        sdf = sdBox(pSquare, vec2(width*f, width));
    else if(shape==2)
        sdf = sdBox(pSquare, vec2(width, width*f));
    else if(shape==3)
        sdf = sdBox(pSquare, vec2(fw)) + delta;
    else if(shape==4)
        sdf = sdRhombus(pSquare, vec2(width*f, width));
    else if(shape==5)
        sdf = sdRhombus(pSquare, vec2(width, width*f));
    else if(shape==6)
        sdf = sdCross(pSquare, vec2(fw, 0.25*fw)) + delta;
    else {
        float ad = 2.0;   // angle divisor, between 2 and n
        int shift;
        if(shape<13){
            ad += f*f*(float(shape)-6.0);
            fw = width;
            delta = 0.0;
            shift = 4;
        }
        else if(shape<19){
            shift = 10;
        }
        else {
            shift = 16;
            ad = float(shape-shift);
        }
        sdf = sdStar(pSquare, fw, shape-shift, ad) + delta;
    }

    vec2 alpha = smoothstep(pixelSize, -pixelSize, sdf + vec2(0.0, outlineWidth));
    outColor = mix(outlineColor, fillColor, alpha.y); // at 0: completely outlineColor, at1: completely fillColor
    outColor.a *= alpha.x;
}