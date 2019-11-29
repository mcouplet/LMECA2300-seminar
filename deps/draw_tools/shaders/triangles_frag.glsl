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

layout (std140) uniform worldBlock
{
	vec2 resolution;
	vec2 translate;
	float zoom;
	// float rotation;
};

layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;


in vec3 bary;
flat in float pixelSize;

out vec4 outColor;


void main()
{
	outColor = fillColor;
	return;
	// vec2 sdf = vec2(0.0, 0.0 - outlineWidth + step(outlineWidth, 0.0))
	//           - min(bary.x, min(bary.y, bary.z));

	// vec2 alpha = smoothstep(0, -2*pixelSize, sdf);
	// outColor = mix(outlineColor, fillColor, alpha.y); // at 0: completely outlineColor, at1: completely fillColor
	// outColor.a *= alpha.x;
}