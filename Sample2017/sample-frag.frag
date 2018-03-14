#version 400
#extension GL_ARB_separate_shader_objects  : enable
#extension GL_ARB_shading_language_420pack : enable


layout( std140, set = 1, binding = 0 ) uniform lightBuf
{
        vec4 uLightPos;
} Light;


layout( std140, set = 2, binding = 0 ) uniform miscBuf
{
	float uTime;
	int   uMode;
} Misc;

// opaque must be outside of a uniform block:
// also, can't specify packing
layout( set = 3, binding = 0 ) uniform sampler2D uSampler;

layout ( location = 0 ) in vec3 vNormal;
layout ( location = 1 ) in vec3 vColor;
layout ( location = 2 ) in vec2 vTexCoord;

layout ( location = 0 ) out vec4 fFragColor;

void
main( )
{
	vec3 rgb;
	switch( Misc.uMode )
	{
		case 0:
			rgb = vColor;
			break;

		case 1:
			rgb = texture( uSampler, vTexCoord ).rgb;
			break;

		default:
			rgb = vec3( 1., 1., 0. );
	}

	fFragColor = vec4( rgb, 1. );
}
