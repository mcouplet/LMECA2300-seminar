#version 150 core

layout (std140) uniform objectBlock
{
    vec4 fillColor;
    vec4 outlineColor;
    vec2 outlineShift; // xy: outline shift x and y
    vec2 localPos;
    vec2 textwidth;
    float boldness; // width here is the boldness
    float outlineWidth;
    // float rotation;
    int space_type; // 0: normal sizes, 1: size in pixels, 2: size in pixels without translation
};

// layout (std140) uniform worldBlock
// {
//     vec2 resolution;
//     vec2 scale;
//     vec2 translate;
//     // float rotation;
// };

uniform sampler2D fontTex;

in vec2 texCoord;
out vec4 outColor;


void main(void)
{
    // vec2 dpdx = dFdx(texCoord);
    // vec2 dpdy = dFdy(texCoord);
    // vec4 texValue = textureGrad(fontTex, texCoord, dpdx, dpdy);
    vec4 texValue = texture(fontTex, texCoord);
    float sdf = texValue.r;
    vec2 normal = texValue.gb - 0.5;
    // vec2 dxdy = vec2(dFdx(sdf), dFdy(sdf));
    float sdfWidth = fwidth(sdf); // length(dxdy);

    // Outline
    float shift = outlineWidth;

    // we add the components of the outlineShift in the direction of the gradient.
    if(outlineShift.xy!=vec2(0.0)) { shift += dot(outlineShift.xy, normal); }

    float glyph_center = 0.5 - 0.25*boldness;
    float outline_center = glyph_center + 0.25*shift;

    // glyph_center = min(glyph_center, outline_center);
    float opacity = smoothstep(glyph_center-sdfWidth, glyph_center+sdfWidth, sdf); // ~2 pixels antialising
    float mu = smoothstep(outline_center-sdfWidth, outline_center+sdfWidth, sdf);
    outColor = mix(outlineColor, fillColor, mu); // at 0: completely outlineColor, at1: completely fillColor
    outColor.a *= opacity;
}