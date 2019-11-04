#version 150 core

layout (std140) uniform objectBlock
{
    vec4 fillColor;
    vec4 outlineColor;
    vec2 otherParam; // xy: outline shift x and y
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
};

uniform sampler2D fontTex;

in vec2 texCoord;
out vec4 outColor;


void main(void)
{
    float sdf = texture(fontTex, texCoord).r;

    // vec2 width_scaling = localScale*(space_type!=0 ? 2.0/resolution : scale);

    float glyph_center = 0.5 - 0.2*width;
    vec2 dxdy = vec2(dFdx( sdf ), dFdy( sdf ));
    // float width = fwidth(sdf);
    float sdfWidth = length(dxdy);
    // float sdfWidth = max(abs(dxdy.x), abs(dxdy.y));
    float opacity = smoothstep(glyph_center-sdfWidth, glyph_center+sdfWidth, sdf); // ~2 pixels antialising

    // Outline
    float shift = 1.0;
    if(otherParam.xy!=vec2(0.0) && sdfWidth>0.00006) { shift += dot(otherParam.xy, dxdy/sdfWidth); }
    float outline_center = glyph_center + 0.2*outlineWidth*shift;
    float mu = smoothstep(outline_center-sdfWidth, outline_center+sdfWidth, sdf);
    outColor = mix(outlineColor, fillColor, mu); // at 0: completely outlineColor, at1: completely fillColor
    outColor.a *= opacity;
}