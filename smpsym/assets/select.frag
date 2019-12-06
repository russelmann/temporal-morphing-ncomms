#version 150

//uniform float uBrightness;

in vec4 color;

out vec4 oFragColor;

void main()
{
	//oFragColor.rgb = mix( vFaceColor, vEdgeColor, fEdgeIntensity ) * uBrightness;
	oFragColor = color;
}
