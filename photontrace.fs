#version 120


uniform sampler2D RandomTexture;

uniform int MaxPathLength;

uniform sampler2D TextureLightSources;
uniform float LightSourceStride;
uniform float LightSummedArea;
uniform sampler3D VolumeTextureTextures;

uniform sampler2D RayDirectionTexture;
uniform sampler2D RayOriginTexture;
uniform sampler2D TexturePolygons;
uniform vec2 PolygonDataStride;

uniform sampler2D TextureBVH;
uniform samplerCube CubeTextureBBoxRootIndices;

uniform sampler2D TextureMaterials;

uniform sampler2D PhotonIntersectionTexture;
uniform sampler2D PhotonPositionTexture;
uniform sampler2D PhotonFluxTexture;
uniform sampler2D PhotonDirectionTexture;

uniform vec4 SceneBSphere;

uniform float MaterialStride;
uniform float MaterialNumRcp;

uniform vec4 OffsetToBBoxMinMax;


const float eps = 1e-6;


struct Ray
{
	vec3 org;
	vec3 dir;
};

struct Intersection
{
	float t;
	vec3 pos;
	vec3 nrm;
	int brdf;
	vec3 col;
	float eta;
	vec3 gnrm;
	vec2 tex;
	float matid;
	float g;
};


vec3 onb(const vec3 x, const vec3 n)
{
	vec3 u, w, v;
	v = n;
	 
	if (n.z < -0.9999999)
	{
		u = vec3(0.0, -1.0, 0.0); 
		w = vec3(-1.0, 0.0, 0.0);
	}
	else
	{
		float a = 1.0 / (1.0 + n.z); 
		float b = -n.x * n.y * a;
		u = vec3(1.0 - n.x * n.x * a, b, -n.x); 
		w = vec3(b, 1.0 - n.y * n.y * a, -n.y);
	}
	return (x.x * u + x.y * v + x.z * w);
}


float GPURnd(inout vec4 n)
{
	// Based on the post http://gpgpu.org/forums/viewtopic.php?t=2591&sid=17051481b9f78fb49fba5b98a5e0f1f3
	// (The page no longer exists as of March 17th, 2015. Please let me know if you see why this code works.)
	const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);
	const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);
	const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);
	const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);

	vec4 beta = floor(n / q);
	vec4 p = a * (n - beta * q) - beta * r;
	beta = (sign(-p) + vec4(1.0)) * vec4(0.5) * m;
	n = (p + beta);

	return fract(dot(n / m, vec4(1.0, -1.0, 1.0, -1.0)));
}


void next(inout vec2 TriangleIndex, const float offset)
{
	TriangleIndex.y = TriangleIndex.y + PolygonDataStride.y * floor(TriangleIndex.x + PolygonDataStride.x * offset);
	TriangleIndex.x = fract(TriangleIndex.x + PolygonDataStride.x * offset);
}


vec4 LastIntersection;
Intersection raytrace(const Ray ray)
{
	Intersection Result;

	Result.t = 1.0e+30;
	Result.nrm = vec3(0.0);
	Result.pos = vec3(0.0);
	Result.col = vec3(0.0);

	vec3 RayDirection = ray.dir;
	vec3 RayOrigin = ray.org;

	vec4 isect = vec4(-1.0, -1.0, -1.0, 1.0e30);
	vec3 Barycentric;

	// RootIndex is the offset to the root
	vec4 RootIndex = textureCube(CubeTextureBBoxRootIndices, RayDirection);
	vec4 BBoxIndex = RootIndex;

	vec3 RayDirectionRcp = vec3(1.0) / RayDirection;
	vec3 t0 = -RayDirectionRcp * RayOrigin;

	while (true)
	{
		vec4 BBoxMinMaxIndex = OffsetToBBoxMinMax + (BBoxIndex.xyxy - RootIndex.xyxy);
		vec4 BBoxMinTriangleX = texture2D(TextureBVH, BBoxMinMaxIndex.xy);
		vec4 BBoxMaxTriangleY = texture2D(TextureBVH, BBoxMinMaxIndex.zw);

		vec4 BBoxNextIndex = texture2D(TextureBVH, BBoxIndex.xy);

		vec3 BBMinInterval = RayDirectionRcp * BBoxMinTriangleX.xyz + t0;
		vec3 BBMaxInterval = RayDirectionRcp * BBoxMaxTriangleY.xyz + t0;
		vec3 a = min(BBMinInterval, BBMaxInterval);
		vec3 b = max(BBMinInterval, BBMaxInterval);
		float tmin = max(max(a.x, a.y), a.z);
		float tmax = min(min(b.x, b.y), b.z);

		bool BBHit = (tmin <= tmax) && (tmin <= isect.w) && (tmax >= 0.0);

		// texture fetches outside seems to be better...
		// read triangle vertices
		vec2 TriangleStartIndex = vec2(BBoxMinTriangleX.w, BBoxMaxTriangleY.w);
		vec2 TriangleIndex = TriangleStartIndex;

		// (px, py, pz, tx)
		// (nx, ny, sgn(nz) * (matid + 1), ty)
		vec3 V0 = texture2D(TexturePolygons, TriangleIndex).xyz; next(TriangleIndex, 1.0); 
		vec3 V1 = texture2D(TexturePolygons, TriangleIndex).xyz; next(TriangleIndex, 1.0); 
		vec3 V2 = texture2D(TexturePolygons, TriangleIndex).xyz; next(TriangleIndex, 1.0); 
		float MaterialIndex = abs(texture2D(TexturePolygons, TriangleIndex).z) - 1.0;
		//int BRDF = BRDFs[int(MaterialIndex)];
		int BRDF =  int(texture2D(TextureMaterials, vec2((MaterialIndex + 0.5 + 0.25) * MaterialStride, 0.0)).x);

		// perform ray-triangle intersection if it is a leaf node
		if ((BBoxMinTriangleX.w >= 0.0) && BBHit)
		{
			// ray triangle intersection
			vec3 p0 = V0;
			vec3 e0 = V1 - V0;
			vec3 e1 = V2 - V0;
			vec3 pv = cross(RayDirection, e1);

			float det = dot(e0, pv);
			{
				vec3 tv = RayOrigin - p0;
				vec3 qv = cross(tv, e0);

				vec4 uvt;
				uvt.x = dot(tv, pv);
				uvt.y = dot(RayDirection, qv);
				uvt.z = dot(e1, qv);
				uvt.xyz = uvt.xyz / det;
				uvt.w = 1.0 - uvt.x - uvt.y;

				if (all(greaterThanEqual(uvt, vec4(0.0))) && (uvt.z < isect.a) && (BRDF != -1))
				{
					Barycentric = uvt.ywx;
					isect = vec4(TriangleStartIndex, MaterialIndex, uvt.z);
				}
			}
		}

		if (BBHit)
		{
			// hit
			BBoxIndex.xy = BBoxNextIndex.xy;
		}
		else
		{
			// miss
			BBoxIndex.xy = BBoxNextIndex.wz;
		}
		if (BBoxIndex.x < 0.0) break;
	};

	if (isect.x >= 0.0)
	{
		vec2 TriangleIndex = isect.xy;

		vec4 V0 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);
		vec4 V1 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);
		vec4 V2 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);

		vec4 N0 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);
		vec4 N1 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);
		vec4 N2 = texture2D(TexturePolygons, TriangleIndex); 

		N0.z = sign(N0.z) * sqrt(abs(1.0 - N0.x * N0.x - N0.y * N0.y)); 
		N1.z = sign(N1.z) * sqrt(abs(1.0 - N1.x * N1.x - N1.y * N1.y)); 
		N2.z = sign(N2.z) * sqrt(abs(1.0 - N2.x * N2.x - N2.y * N2.y)); 

		Result.t = isect.w;
		Result.pos = V2.xyz * Barycentric.x + V0.xyz * Barycentric.y + V1.xyz * Barycentric.z;
		Result.nrm = normalize(N2.xyz * Barycentric.x + N0.xyz * Barycentric.y + N1.xyz * Barycentric.z);
		Result.gnrm = normalize(cross(V1.xyz - V0.xyz, V2.xyz - V0.xyz));

		vec2 T0 = vec2(V0.w, N0.w);
		vec2 T1 = vec2(V1.w, N1.w);
		vec2 T2 = vec2(V2.w, N2.w);
		Result.tex = T2 * Barycentric.x + T0 * Barycentric.y + T1 * Barycentric.z;

		Result.col = texture2D(TextureMaterials, vec2((isect.z + 0.0 + 0.25) * MaterialStride, 0.0)).xyz;
		Result.brdf = int(texture2D(TextureMaterials, vec2((isect.z + 0.5 + 0.25) * MaterialStride, 0.0)).x);
		Result.g = texture2D(TextureMaterials, vec2((isect.z + 0.5 + 0.25) * MaterialStride, 0.0)).y;
		Result.eta = texture2D(TextureMaterials, vec2((isect.z + 0.0 + 0.25) * MaterialStride, 0.0)).w;
		Result.matid = isect.z;
	}

	LastIntersection = isect;
	return Result;
}


vec3 glossy_reflect(const vec3 d, const vec3 n, const float g, inout vec4 rndv)
{
	float a = 2.0 / (g + 1.0);
	vec3 r = normalize((1.0 - a) * reflect(d, n) + a * n);

	float rnd1 = GPURnd(rndv);
	float rnd2 = GPURnd(rndv);

	float temp1 = 2.0 * 3.141592 * rnd1;
	float temp2 = sqrt(1.0 - pow(rnd2, 2.0 / (g + 1.0)));
	vec3 v = vec3(sin(temp1) * temp2, pow(rnd2, 1.0 / (g + 1.0)), cos(temp1) * temp2);

	vec3 result = normalize(onb(v, r));
	if (dot(result, n) < 0.0)
	{
		result = reflect(result, n);
	}

	return result;
}


void GenerateIBLSample(inout vec4 rndv, out vec3 flux, out vec3 dir, out vec3 org)
{
	float SceneProjectedArea = 4 * 3.141592 * SceneBSphere.w * SceneBSphere.w;

	float temp1 = 2.0 * acos(sqrt(1.0 - GPURnd(rndv)));
	float temp2 = 2.0 * 3.141592 * GPURnd(rndv);
	dir.x = sin(temp1) * cos(temp2);
	dir.y = cos(temp1);
	dir.z = sin(temp1) * sin(temp2);

	if (GPURnd(rndv) > 0.5)
	{
		dir = vec3(1.0);
		dir.y = 2.0 * dir.y;
		dir.z = -dir.z;
		dir = normalize(dir);
		flux = vec3(1.0, 1.0, 0.7) * SceneProjectedArea;
	}
	else
	{
		flux = vec3(0.7, 0.7, 1.0) * SceneProjectedArea;
	}

	float radius = sqrt(GPURnd(rndv));
	float theta = 2.0 * 3.141592 * GPURnd(rndv);
	vec3 temp = vec3(radius * cos(theta), 0.0, radius * sin(theta));
	org = onb(temp * SceneBSphere.w, dir) + dir * SceneBSphere.w + SceneBSphere.xyz;
	dir = -dir;
}



void GenerateLightSourceSample(inout vec4 rndv, out vec3 flux, out vec3 dir, out vec3 org)
{
	// generate uniform random values
	float Rnd = GPURnd(rndv);
	float rndv1 = Rnd;

	// binary search
	int low = 0;
	int high = int(1.0 / LightSourceStride + 0.5);
	vec2 LightIndex;
	while (low < high) 
	{
		int mid = low + ((high - low) / int(2));  
		vec2 midv = vec2((float(mid) - 0.5) * LightSourceStride, 0.0);
		vec3 tmpv = texture2D(TextureLightSources, midv).xyz;
		LightIndex = tmpv.xy;

		if (tmpv.z < Rnd)
		{
			low = mid + 1; 
		}
		else
		{
			high = mid; 
		}
	}
	{
		vec2 midv = vec2((float(high) - 0.5) * LightSourceStride, 0.0);
		vec3 tmpv = texture2D(TextureLightSources, midv).xyz;
		LightIndex = tmpv.xy;
	}

	// calculate paramters
	float rndv2 = GPURnd(rndv);
	float t = sqrt(1.0 - rndv2);
	float s = GPURnd(rndv);
	float rndv3 = s;

	float a = 1.0 - t;
	float b = (1.0 - s) * t;
	float c = s * t;

	// interpolate the position and the normal
	vec2 TriangleIndex = LightIndex;

	vec4 V0 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);
	vec4 V1 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);
	vec4 V2 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);

	vec4 N0 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);
	vec4 N1 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);
	vec4 N2 = texture2D(TexturePolygons, TriangleIndex); 
	float MaterialIndex = abs(N0.z) - 1.0;

	N0.z = sign(N0.z) * sqrt(abs(1.0 - N0.x * N0.x - N0.y * N0.y)); 
	N1.z = sign(N1.z) * sqrt(abs(1.0 - N1.x * N1.x - N1.y * N1.y)); 
	N2.z = sign(N2.z) * sqrt(abs(1.0 - N2.x * N2.x - N2.y * N2.y)); 

	flux = texture2D(TextureMaterials, vec2((MaterialIndex + 0.0 + 0.25) * MaterialStride, 0.0)).xyz;

	flux *= vec3(LightSummedArea * 3.141592);

	vec2 T0 = vec2(V0.w, N0.w);
	vec2 T1 = vec2(V1.w, N1.w);
	vec2 T2 = vec2(V2.w, N2.w);
	vec2 T = T2 * c + T0 * a + T1 * b;

	if (abs(T.x) < 1e+10) 
	{
		flux *= texture3D(VolumeTextureTextures, vec3(T, MaterialNumRcp * (MaterialIndex + 0.5 + 0.25))).rgb;
	}

	org = V0.xyz * a + V1.xyz * b + V2.xyz * c;

	// uniform light
	vec3 nrm = normalize(N0.xyz * a + N1.xyz * b + N2.xyz * c);
	vec3 gnrm = normalize(cross(V1.xyz - V0.xyz, V2.xyz - V0.xyz));

	org = org + gnrm * 1.0e-5;

	vec2 rnd;
	rnd.x = GPURnd(rndv);
	rnd.y = GPURnd(rndv);
	rnd.x = 2.0 * 3.141592 * rnd.x;
	rnd.y = sqrt(rnd.y);
	dir = onb(vec3(sin(rnd.x) * rnd.y, sqrt(1.0 - rnd.y * rnd.y), cos(rnd.x) * rnd.y), nrm);
}


float Fresnel(in vec3 incom, in vec3 normal, in float index_internal, in float index_external) 
{
	float eta = index_internal / index_external;
	float cos_theta1 = dot(incom, normal);
	float cos_theta2 = 1.0 - (eta * eta) * (1.0 - cos_theta1 * cos_theta1);

	if (cos_theta2 < 0.0)
	{
		return 1.0;
	}
	else
	{
		cos_theta2 = sqrt(cos_theta2);
		float fresnel_rs = (index_internal * cos_theta1 - index_external * cos_theta2) / (index_internal * cos_theta1 + index_external * cos_theta2);
		float fresnel_rp = (index_internal * cos_theta2 - index_external * cos_theta1) / (index_internal * cos_theta2 + index_external * cos_theta1);
		return (fresnel_rs * fresnel_rs + fresnel_rp * fresnel_rp) * 0.5;
	}
}


void main()
{
	vec2 PixelIndex = gl_TexCoord[0].st; 

	// state of the random number generator
	vec4 rndv = texture2D(RandomTexture, PixelIndex);

	// read previous intersection
	vec4 PhotonPosition = texture2D(PhotonPositionTexture, PixelIndex);
	vec4 PhotonFlux = texture2D(PhotonFluxTexture, PixelIndex);
	vec4 PhotonDirection = texture2D(PhotonDirectionTexture, PixelIndex);
	vec3 PhotonNrm = vec3(PhotonPosition.w, PhotonFlux.w, PhotonDirection.w);

	vec3 flux = abs(PhotonFlux.rgb);
	Ray r;


	bool ContinueTrace;
	float EmittedFlag;

	vec4 PhotonIntersection = texture2D(PhotonIntersectionTexture, PixelIndex);
	PhotonIntersection.z += 1.0;

	if ((PhotonIntersection.x < 0.0) || (PhotonIntersection.z >= MaxPathLength))
	{
		// the last photon trace was terminated, generate a new photon 
		if (LightSummedArea != 0.0)
		{
			GenerateLightSourceSample(rndv, flux, r.dir, r.org);
		}
		else
		{
			GenerateIBLSample(rndv, flux, r.dir, r.org);
		}

		EmittedFlag = 1.0;

		// reset the trace level
		PhotonIntersection.z = 0.0;
		ContinueTrace = true;
	}
	else
	{
		// extra bounce
		Intersection i;
		vec2 TriangleIndex = PhotonIntersection.xy;

		vec4 V0 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);
		vec4 V1 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);
		vec4 V2 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);

		vec4 N0 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);
		vec4 N1 = texture2D(TexturePolygons, TriangleIndex); next(TriangleIndex, 1.0);
		vec4 N2 = texture2D(TexturePolygons, TriangleIndex); 

		i.gnrm = normalize(cross(V1.xyz - V0.xyz, V2.xyz - V0.xyz));


		vec4 Barycentric;
		vec3 t0, t1, t2;

		t0 = V1.xyz - V1.xyz;
		t1 = PhotonPosition.xyz - V0.xyz;
		t2 = cross(t0, t1);
		Barycentric.x = length(t2);

		t0 = V2.xyz - V1.xyz;
		t1 = PhotonPosition.xyz - V1.xyz;
		t2 = cross(t0, t1);
		Barycentric.y = length(t2);

		t0 = V0.xyz - V2.xyz;
		t1 = PhotonPosition.xyz - V2.xyz;
		t2 = cross(t0, t1);
		Barycentric.z = length(t2);

		Barycentric.w = Barycentric.x + Barycentric.y + Barycentric.z;
		Barycentric = Barycentric / Barycentric.w;

		vec2 T0 = vec2(V0.w, N0.w);
		vec2 T1 = vec2(V1.w, N1.w);
		vec2 T2 = vec2(V2.w, N2.w);
		i.tex = T2 * Barycentric.x + T0 * Barycentric.y + T1 * Barycentric.z;


		// last intersection
		i.t = PhotonIntersection.w; // later it will be used for distance based attenuation
		i.pos = PhotonPosition.xyz;
		i.nrm = PhotonNrm;
		float MaterialIndex = abs(N0.z) - 1.0;
		i.col = texture2D(TextureMaterials, vec2((MaterialIndex + 0.25) * MaterialStride, 0.0)).xyz;
		//i.brdf = BRDFs[int(MaterialIndex + 0.25)];
		i.brdf = int(texture2D(TextureMaterials, vec2((MaterialIndex + 0.5 + 0.25) * MaterialStride, 0.0)).x);
		i.eta = texture2D(TextureMaterials, vec2((MaterialIndex + 0.25) * MaterialStride, 0.0)).w; 
		i.g = texture2D(TextureMaterials, vec2((MaterialIndex + 0.5 + 0.25) * MaterialStride, 0.0)).y; 

		if (abs(i.tex.x) < 1e+10) 
		{
			i.col *= texture3D(VolumeTextureTextures, vec3(i.tex, MaterialNumRcp * (MaterialIndex + 0.5 + 0.25))).rgb;
		}

		if (i.brdf == 0)
		{
			// matte
			float r0 = 2.0 * 3.141592 * GPURnd(rndv);
			float r1 = sqrt(GPURnd(rndv));
			vec3 v = vec3(sin(r0) * r1, sqrt(1.0 - r1 * r1), cos(r0) * r1);

			r.org = i.pos + eps * i.gnrm;
			r.dir = onb(v, i.nrm);
			if (dot(r.dir, i.gnrm) < 0.0) r.dir = -r.dir;
		}
		else if ((i.brdf == 1) || (i.brdf == 4))
		{
			// metal
			r.org = i.pos + eps * i.gnrm;
			if (i.brdf == 4) i.nrm = glossy_reflect(i.nrm, i.nrm, 1.0 / pow((1.0 - i.g) * 0.5, 2.71828), rndv);
			r.dir = reflect(PhotonDirection.xyz, i.nrm);

			if (dot(r.dir, i.gnrm) < 0.0) r.dir = -r.dir;
		}
		else if ((i.brdf == 2) || (i.brdf == 5))
		{
			if (i.brdf == 5) i.nrm = glossy_reflect(i.nrm, i.nrm, 1.0 / pow((1.0 - i.g) * 0.5, 2.71828), rndv);

			// dielectric
			float ln = dot(i.nrm, PhotonDirection.xyz);
			float eta = i.eta;

			if (ln < 0.0)
			{
				// ray is going in
				float Re = Fresnel(-PhotonDirection.xyz, i.nrm, 1.0, eta);
				if (GPURnd(rndv) < Re)
				{
					// specular reflection
					r.org = i.pos + eps * i.gnrm;
					r.dir = reflect(PhotonDirection.xyz, i.nrm);
					i.col = vec3(1.0);
					if (dot(r.dir, i.gnrm) < 0.0) r.dir = -r.dir;
				}
				else
				{
					// specular refraction
					r.org = i.pos - eps * i.gnrm;
					r.dir = refract(PhotonDirection.xyz, i.nrm, 1.0 / eta);
					if (dot(r.dir, i.gnrm) > 0.0) r.dir = -r.dir;
				}
			}
			else
			{
				float Re = Fresnel(-PhotonDirection.xyz, -i.nrm, eta, 1.0);
				if (GPURnd(rndv) < Re)
				{
					// specular reflection
					r.org = i.pos - eps * i.gnrm;
					r.dir = reflect(PhotonDirection.xyz, -i.nrm);
					if (dot(r.dir, -i.gnrm) < 0.0) r.dir = -r.dir;
				}
				else
				{
					// specular refraction
					r.org = i.pos + eps * i.gnrm;
					r.dir = refract(PhotonDirection.xyz, -i.nrm, eta);
					if (dot(r.dir, -i.gnrm) > 0.0) r.dir = -r.dir;
				}
			}
		}
		else if ((i.brdf == 3) || (i.brdf == 6))
		{
			if (i.brdf == 6) i.nrm = glossy_reflect(i.nrm, i.nrm, 1.0 / pow((1.0 - i.g) * 0.5, 2.71828), rndv);

			// plastic
			float ln = -abs(dot(i.nrm, PhotonDirection.xyz));
			float eta = i.eta;
			float Re = Fresnel(-PhotonDirection.xyz, i.nrm, 1.0, eta);
			if (GPURnd(rndv) < Re)
			{
				// specular reflection (assume that the color of the coating is 1.0)
				r.org = i.pos + eps * i.gnrm;
				r.dir = reflect(PhotonDirection.xyz, i.nrm);
				i.col = vec3(1.0);
				if (dot(r.dir, i.gnrm) < 0.0) r.dir = -r.dir;
			}
			else
			{
				// matte
				float r0 = 2.0 * 3.141592 * GPURnd(rndv);
				float r1 = sqrt(GPURnd(rndv));
				vec3 v = vec3(sin(r0) * r1, sqrt(1.0 - r1 * r1), cos(r0) * r1);

				r.org = i.pos + eps * i.gnrm;
				r.dir = onb(v, i.nrm);
				if (dot(r.dir, i.gnrm) < 0.0) r.dir = -r.dir;
			}
		}

		// Russian roulette
		float p = max(max(i.col.r, i.col.g), i.col.b);
		if (p < GPURnd(rndv))
		{
			// photon is terminated
			ContinueTrace = false;
		}
		else
		{
			// continue tracing
			ContinueTrace = true;
			i.col = i.col / p;
			flux = flux * i.col;
		}
		EmittedFlag = 0.0;
	}



	if (ContinueTrace)
	{
		Intersection i = raytrace(r);

		if ((i.t == 1.0e+30) || ((i.brdf != 0) && (i.brdf != 6) && (i.brdf != 3)))
		{
			// no intersection, invalidate the photon
			flux = -flux;
		}

		float gfactor = min(abs(dot(r.dir, i.nrm) / dot(r.dir, i.gnrm)), sqrt(5.0));
		flux *= gfactor;

		LastIntersection.z = PhotonIntersection.z;
		gl_FragData[0] = vec4(i.pos, i.nrm.x);
		gl_FragData[1] = vec4( flux, i.nrm.y);
		gl_FragData[2] = vec4(r.dir, i.nrm.z);
		gl_FragData[3] = rndv;
		gl_FragData[4] = LastIntersection;
		gl_FragData[5] = vec4(EmittedFlag);
	}
	else
	{
		gl_FragData[0] = vec4(0.0);
		gl_FragData[1] = vec4(vec3(-1.0), 0.0);
		gl_FragData[2] = vec4(0.0);
		gl_FragData[3] = rndv;
		gl_FragData[4] = vec4(-1.0, -1.0, 0.0, 1e+30);
		gl_FragData[5] = vec4(EmittedFlag);
	}
}
