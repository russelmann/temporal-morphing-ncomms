#version 150

uniform mat4	ciModelViewProjection;
in vec4		ciPosition;
in vec4		ciColor;
//in vec2		ciTexCoord0;

out vec4 color;

void main()
{
	color = ciColor;
	gl_Position = ciModelViewProjection * ciPosition;
}
