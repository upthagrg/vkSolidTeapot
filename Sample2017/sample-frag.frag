#version 400
#extension GL_ARB_separate_shader_objects  : enable
#extension GL_ARB_shading_language_420pack : enable


layout( std140, set = 1, binding = 0 ) uniform lightBuf
{
	float uKa;
	float uKd;
	float uKs;
	vec4 uLightPos;
	vec3 uLightSpecularColor;
	float uShininess;
} Light;


layout( std140, set = 2, binding = 0 ) uniform miscBuf
{
	float uTime;
	int   uMode;
	int   uLit;
	int   uModulate;
} Misc;

// opaque must be outside of a uniform block:
// also, can't specify packing
layout( set = 3, binding = 0 ) uniform sampler2D uSampler;

layout ( location = 0 ) in vec3 vColor;
layout ( location = 1 ) in vec2 vTexCoord;
layout ( location = 2 ) in vec3 vN;
layout ( location = 3 ) in vec3 vL;
layout ( location = 4 ) in vec3 vE;

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
//if using lighting
//if( ????? != 0 ) 
	if( Misc.uLit != 0 )  
	{
		vec3 normal = normalize(vN);
		vec3 light  = normalize(vL);
		vec3 eye    = normalize(vE);

		//vec3 ambient = ?????
		vec3 ambient;

		if(Misc.uMode == 1){//textured
			if(Misc.uModulate == 0){
				ambient = texture( uSampler, vTexCoord ).rgb;
			}
			else{
				ambient = vColor.rgb;
			}
		}
		else{
			ambient = vColor.rgb;
		}
		float d = 0.;
		if( dot(normal,light) > 0. )
		{
			//d = ?????
			d = 1.;
		}
		//vec3 diffuse = ?????
		vec3 diffuse;
		if(Misc.uMode==1){//textured
			diffuse = d * max(dot(light,normal),0.)* texture( uSampler, vTexCoord ).rgb;
		}
		else{
			diffuse = d * max(dot(light,normal),0.)*vColor.rgb;
		}
		float s = 0.;
		if( dot(normal,light) > 0. )		// only do specular if the light can see the point
		{
			vec3 ref = reflect( -eye, normal );
			//s = ?????
			s = 1.;
		}
		//vec3 specular = ?????
		vec3 specular = Light.uLightSpecularColor.rgb * pow(max(dot(reflect(-light,normal),eye),0.), Light.uShininess) * s;

		rgb = (Light.uKa*ambient) + (Light.uKd*diffuse) + (Light.uKs*specular);
	}

	fFragColor = vec4( rgb, 1. );
}