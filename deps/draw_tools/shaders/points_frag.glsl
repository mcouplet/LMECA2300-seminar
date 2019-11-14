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

// p is the coordinate, s is the size
float sdCross(vec2 p, vec2 s, float z)
{
    vec2 zs = z*s;
    p = abs(p); p = (p.y>p.x) ? p.yx : p.xy;
    vec2  q = p - zs;
    float k = max(q.y,q.x);
    vec2  w = (k>0.0) ? q : vec2(zs.y-p.x,-k);
    return sign(k)*length(max(w,0.0)) - s.x*(1.0-z);
}

float sdCircle(vec2 p, float r) {
    return length(p) - r;
}


float sdPentagon(vec2 p, float r, float z)
{
    const vec3 k = vec3(0.809016994,0.587785252,0.726542528);
    float zr = k.x*z*r;
    p.x = abs(p.x);
    p -= 2.0*min(dot(vec2(-k.x,k.y),p),0.0)*vec2(-k.x,k.y);
    p -= 2.0*min(dot(vec2( k.x,k.y),p),0.0)*vec2( k.x,k.y);
    p -= vec2(clamp(p.x,-zr*k.z, zr*k.z), zr);
    return length(p)*sign(p.y) - r*(1.0-z);
}

float sdBox(vec2 p, float l, float z)
{

    vec2 d = abs(p);
    d -= vec2(l, l*min(z, 1./z));
    return length(max(d,vec2(0))) + min(max(d.x,d.y),0.0);
}

float ndot(vec2 a, vec2 b ) { return a.x*b.x - a.y*b.y; }

float sdRhombus(vec2 p, vec2 b) 
{
    vec2 q = abs(p);

    float h = clamp( (-2.0*ndot(q,b) + ndot(b,b) )/dot(b,b), -1.0, 1.0 );
    float d = length( q - 0.5*b*vec2(1.0-h,1.0+h) );
    d *= sign( q.x*b.y + q.y*b.x - b.x*b.y );
    
    return d;
}


// signed distance to a centered line
float sdCenteredLine(vec2 p, vec2 b, float l2)
{
    float h = clamp( dot(p,b)/l2, 0.0, 1.0 );
    return length( p - b*h );
}


float sdTripod(vec2 p, float r) {
    p.x = abs(p.x);
    float r2 = r*r;
    float d0 = sdCenteredLine( p, r*vec2( 0.5*sqrt(3.), -0.5), r2 );
    float d1 = sdCenteredLine( p, vec2( 0.0,  r), r2 );
    return min(d0, d1);
}

float sdQuadPod(vec2 p, float r) {
    // p = abs(vec2(p.x+p.y, p.x-p.y));// gives a vertical cross
    p = abs(p); // diagonal cross
    return sdCenteredLine( p, r*vec2( 1., 1.), r*r*2.0 );
}


// much more interesting: allows to draw stars and hexagon and circle !
float sdHexagram(vec2 p, float r)
{
    const vec4 k = vec4(-0.5,0.86602540378,0.57735026919,1.73205080757);
    
    p = abs(p);
    p -= 2.0*min(dot(k.xy,p),0.0)*k.xy;
    p -= 2.0*min(dot(k.yx,p),0.0)*k.yx;
    p -= vec2(clamp(p.x,r*k.z,r*k.w),r);
    return length(p)*sign(p.y);
}



void main( void ) {
    // float z = 1.5*fract(marker*0.1);
    // vec2 sdf = sdCross(pSquare, vec2(width, 0.25*width), z) + vec2(0.0, outlineWidth);

    // float z = 1.0-1.0/(1.0+marker);
    // vec2 sdf = sdPentagon(pSquare, width, z) + vec2(0.0, outlineWidth);

    float z = 2.0*fract(marker*0.5);
    if(z>1.0)
        z = 1.0/(2.0-z);
    vec2 sdf = sdBox(pSquare, width, z) + vec2(0.0, outlineWidth);

    vec2 alpha = smoothstep(pixelSize, -pixelSize, sdf);
    outColor = mix(outlineColor, fillColor, alpha.y); // at 0: completely outlineColor, at1: completely fillColor
    outColor.a *= alpha.x;
}