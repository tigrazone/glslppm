#version 120


uniform sampler2D PhotonPositionTexture;
uniform sampler2D PhotonFluxTexture;

uniform int HashNum;
uniform float GridScale; 
uniform vec3 BBoxMin;

uniform vec4 BufInfo;


vec2 convert1Dto2D(const float t)
{
	return vec2(mod(t, BufInfo.x) + 0.5,  floor(t * BufInfo.z) + 0.5);
}


float hash(const vec3 idx)
{
	// use the same procedure as GPURnd
	vec4 n = vec4(idx, GridScale * 0.5) * 4194304.0 / GridScale;

	const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);
	const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);
	const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);
	const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);

	vec4 beta = floor(n / q);
	vec4 p = a * (n - beta * q) - beta * r;
	beta = (sign(-p) + vec4(1.0)) * vec4(0.5) * m;
	n = (p + beta);

	return floor(fract(dot(n / m, vec4(1.0, -1.0, 1.0, -1.0))) * HashNum);
}


void main()
{
	vec2 PixelIndex = gl_TexCoord[0].st;
	vec3 PhotonPosition = texture2D(PhotonPositionTexture, PixelIndex).xyz;
	vec3 HashIndex = floor((PhotonPosition - BBoxMin) * GridScale);

	// vec4(original 2d index, hashed 2d index)
	vec4 PhotonIndex = vec4(PixelIndex.xy, convert1Dto2D(hash(HashIndex)));

	// ignore invalid photons (move them outside the buffer)
	float PhotonFlux = texture2D(PhotonFluxTexture, PixelIndex).x;
	if (PhotonFlux < 0.0)
	{
		PhotonIndex.zw = vec2(-1.0);
	}

	gl_FragColor = PhotonIndex;
}
