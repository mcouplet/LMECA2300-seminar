#version 150 core

layout (std140) uniform objectBlock
{
    vec4 fillColor;
    vec4 outlineColor;
    vec2 otherParam;
    vec2 localPos;
    vec2 localScale;
    float width;
    float outlineWidth;
    // float rotation;
    int space_type; // 0: normal sizes, 1: size in pixels, 2: size in pixels without translation
};

layout (std140) uniform worldBlock
{
    vec2 resolution;
    vec2 scale;
    vec2 translate;
    // float rotation;
    // float wtime;
};

layout(lines_adjacency) in;
layout(triangle_strip, max_vertices = 4) out;

flat out vec2 p1_screen;
flat out vec2 p2_screen;

vec2 perpendicular(vec2 v) {
    return vec2(-v.y, v.x);
}

void main() {
    vec2 p0=gl_in[0].gl_Position.xy;
    vec2 p1=gl_in[1].gl_Position.xy;
    vec2 p2=gl_in[2].gl_Position.xy;
    vec2 p3=gl_in[3].gl_Position.xy;

    bool first = p0==p1;
    bool last = p2==p3;

    vec2 v1 = normalize(p2 - p1);          // direction of this segment
    vec2 v0 = first?v1:normalize(p1 - p0); // direction of the prev segment
    vec2 v2 = last?v1:normalize(p3 - p2);  // direction of next segment

    // determine the normal of each of the 3 segments (previous, current, next)
    vec2 n0 = perpendicular(v0);
    vec2 n1 = perpendicular(v1);
    vec2 n2 = perpendicular(v2);

    // determine miter lines by averaging the normals of the 2 segments
    vec2 miter0 = normalize(n0 + n1);    // miter at start of current segment
    vec2 miter1 = normalize(n1 + n2);    // miter at end of current segment

    // determine the length of the miter by projecting it onto normal and then inverse it
    float miter0_size = width / dot(miter0, n1);
    float miter1_size = width / dot(miter1, n1);

    // we apply the transformations to get screen coordinates (we worked in a normalized world space)
    vec2 center_scaling = localScale*(space_type==2 ? 2.0/resolution : scale); // 2/resolution: the size of one screen pixel
    vec2 width_scaling = localScale*(space_type!=0 ? 2.0/resolution : scale);
    vec2 center_translate = space_type==2 ? localPos : scale*(localPos + translate);

    p1_screen = center_scaling*p1 + center_translate;
    p2_screen = center_scaling*p2 + center_translate;
    vec2 miter0_screen = miter0_size*width_scaling*miter0;
    vec2 miter1_screen = miter1_size*width_scaling*miter1;
    vec2 q0 = p1_screen + miter0_screen - (first?width_scaling*v0:vec2(0.0));
    vec2 q1 = p1_screen - miter0_screen - (first?width_scaling*v0:vec2(0.0));
    vec2 q2 = p2_screen + miter1_screen + (last?width_scaling*v2:vec2(0.0));
    vec2 q3 = p2_screen - miter1_screen + (last?width_scaling*v2:vec2(0.0));

    vec2 miter1_perpend = perpendicular(miter1_screen);

    // orientation of triangle q0-q2-q3
    float ori = dot(q0 - q3, miter1_perpend);

    // the intersection of (q0, q1) and (q2, q3):   coordinates are given by q0 + inter*(q1-q0)
    float inter = 0.5*ori/dot(miter0_screen, miter1_perpend);

    gl_Position.zw = vec2(0.0, 1.0);

    // avoid crossing lines
    if(inter>0 && inter<1) {
        gl_Position.xy = inter*(q1-q0)+q0;
        EmitVertex();

        if(ori<0) {
            // the upper part is inverted
            gl_Position.xy = q1;
            EmitVertex();

            gl_Position.xy = q3;
            EmitVertex();
        }
        else {
            // the lower part is inverted
            gl_Position.xy = q2;
            EmitVertex();

            gl_Position.xy = q0;
            EmitVertex();
        }
    }
    else {
        gl_Position.xy = q0;
        EmitVertex();

        gl_Position.xy = q1;
        EmitVertex();

        gl_Position.xy = q2;
        EmitVertex();

        gl_Position.xy = q3;
        EmitVertex();
    }

    EndPrimitive();
}