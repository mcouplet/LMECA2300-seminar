#version 150 core

layout (std140) uniform objectBlock
{
    vec4 fillColor;
    vec4 outlineColor;
    vec2 localPos;
    vec2 outlineShift;
    float textheight;
    float boldness;
    float outlineWidth;
    // float rotation;
    int space_type; // 0: normal sizes, 1: size in pixels, 2: size in pixels without translation
};

layout (std140) uniform worldBlock
{
    vec2 resolution;
    vec2 translate;
    float zoom;
    // float rotation;
};

// position of this vertex compared to
// the bottom-left of the left-most letter
in vec2 pos;

// position of this vertex in the font atlas texture
in vec2 tex;

// just forward text to the fragment shader
out vec2 texCoord;

void main()
{
    texCoord = tex;

    vec2 resRatio = min(resolution.x, resolution.y)/resolution;

    // note: conditional based on uniform should not slow the shader down
    vec2 height = pos * textheight;
    if(space_type==0) {
        // classical case
        gl_Position = vec4(resRatio*zoom*(localPos + translate + height), 0.0, 1.0);
    }
    else if(space_type==1) {
        // no scaling of the height is applied
        gl_Position = vec4(resRatio*(zoom*(localPos + translate) + height), 0.0, 1.0);
    }
    else {
        // everything is given in pixel, from the bottom left corner :-)
        // we round the position to the bottom left of pixel
        // as the rasterization is done with centers, we will have
        // beautifully aligned pixels if the user choose a font size
        // that is a multiple of the font size :-)
        // vec2 pixelPos = floor(localPos + height)+0.5;
        vec2 pixelPos = localPos + height;
        gl_Position = vec4(2.0*pixelPos/resolution - 1.0, 0.0, 1.0);
    }
}