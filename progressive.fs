#version 120


uniform sampler2D QueryPositionTexture;
uniform sampler2D QueryEmissionPhotonCountTexture;
uniform sampler2D QueryFluxRadiusTexture;
uniform sampler2D QueryReflectanceTexture;
uniform sampler2D QueryNormalTexture;
uniform sampler2D QueryIntersectionTexture;

uniform sampler2D HashedPhotonTexture;
uniform sampler2D PhotonFluxTexture;
uniform sampler2D PhotonPositionTexture;
uniform sampler2D PhotonDirectionTexture;
uniform sampler2D PhotonCorrectionTexture;

uniform vec4 BufInfo;

uniform float HashNum;
uniform float GridScale;
uniform vec3 BBoxMin;
uniform vec3 HashMax;

uniform float Alpha;


float hash(const vec3 idx)
{
	// use the same procedure as GPURnd
	// it is the same as the one in hash.fs
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


vec2 convert1Dto2D(const float t)
{
	return vec2(mod(t, BufInfo.x) + 0.5, floor(t * BufInfo.z) + 0.5);
}


vec3 Flux = vec3(0.0);
float PhotonCount = 0.0;
void AccumulatePhotons(const vec3 QueryPosition, const vec3 QueryNormal, const float QueryRadius, const vec3 HashIndex)
{
	// get the photon from the hash buffer
	vec2 HashedPhotonIndex = convert1Dto2D(hash(HashIndex)) * BufInfo.zw;
	vec2 PhotonIndex = texture2D(HashedPhotonTexture, HashedPhotonIndex).xy;

	vec4 PhotonFlux = texture2D(PhotonFluxTexture, PhotonIndex);
	vec4 PhotonPosition = texture2D(PhotonPositionTexture, PhotonIndex);
	vec4 PhotonDirection = texture2D(PhotonDirectionTexture, PhotonIndex);
	vec3 PhotonNormal = vec3(PhotonPosition.w, PhotonFlux.w, PhotonDirection.w);

	// make sure that the photon is actually in the given grid cell
	vec3 RangeMin = HashIndex / GridScale + BBoxMin;
	vec3 RangeMax = (HashIndex + vec3(1.0)) / GridScale + BBoxMin;
	if (all(greaterThan(PhotonPosition.xyz, RangeMin)) && all(lessThan(PhotonPosition.xyz, RangeMax)))
	{
		// photon projection as in "Diffusion-Based Photon Mapping" by L. Schjoeth.
		vec3 PositionDifference = PhotonPosition.xyz - QueryPosition;
		float t = dot(QueryNormal, PositionDifference) / dot(QueryNormal, PhotonDirection.xyz);
		PhotonPosition.xyz = PhotonPosition.xyz + t * PhotonDirection.xyz;

		if ((length(PositionDifference) < QueryRadius) && (dot(QueryNormal, -PhotonDirection.xyz) > 0.0)) 
		{
			float PhotonCorrection = texture2D(PhotonCorrectionTexture, HashedPhotonIndex).x;
			Flux += PhotonFlux.rgb * PhotonCorrection;
			PhotonCount += PhotonCorrection;
		}
	}
}


void main()
{
	vec2 PixelIndex = gl_TexCoord[0].st; 

	vec4 QueryPosition = texture2D(QueryPositionTexture, PixelIndex);
	vec4 QueryNormal = texture2D(QueryNormalTexture, PixelIndex);
	vec4 QueryFluxRadius = texture2D(QueryFluxRadiusTexture, PixelIndex);
	vec4 QueryEmissionPhotonCount = texture2D(QueryEmissionPhotonCountTexture, PixelIndex);
	vec3 QueryFlux = QueryFluxRadius.xyz;
	vec4 QueryReflectance = texture2D(QueryReflectanceTexture, PixelIndex);
	float QueryPhotonCount = QueryEmissionPhotonCount.w;
	float QueryRadius = QueryFluxRadius.w;

	vec3 RangeMin = max((QueryPosition.xyz - vec3(QueryRadius) - BBoxMin) * GridScale, vec3(0.0));
	vec3 RangeMax = min((QueryPosition.xyz + vec3(QueryRadius) - BBoxMin) * GridScale, HashMax);

	vec4 QueryIntersection = texture2D(QueryIntersectionTexture, PixelIndex);
	if (QueryIntersection.x >= 0.0)
	{
		vec3 ii = floor(abs(QueryPosition.xyz - BBoxMin) * GridScale - 0.5); 
		for (int iz = 0; iz <= 1; iz++)
		{
			for (int iy = 0; iy <= 1; iy++)
			{
				for (int ix = 0; ix <= 1; ix++)
				{
					AccumulatePhotons(QueryPosition.xyz, QueryNormal.xyz, QueryRadius, ii + vec3(ix, iy, iz));
				}
			}
		}

		// BRDF (assumes that we stop at Lambertian - we should use the BRDF there in general.)
		Flux *= (QueryReflectance.rgb / 3.141592);

		// progressive density estimation
		float g = min((QueryPhotonCount + Alpha * PhotonCount) / (QueryPhotonCount + PhotonCount), 1.0);
		QueryRadius = QueryRadius * sqrt(g);
		QueryPhotonCount = QueryPhotonCount + PhotonCount * Alpha;
		QueryFlux = (QueryFlux + Flux) * g;
	}

	gl_FragData[0] = vec4(QueryFlux, QueryRadius);
	gl_FragData[1] = vec4(QueryEmissionPhotonCount.rgb, QueryPhotonCount);
}
