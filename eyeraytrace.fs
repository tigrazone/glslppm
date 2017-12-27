#version 120


uniform sampler2D QueryEmissionPhotonCountTexture;
uniform sampler2D RandomTexture;
uniform sampler3D VolumeTextureTextures;

uniform float FocalLength;
uniform int MaxPathLength;
uniform float ApertureSize;
uniform float LightSummedArea;

uniform int NumEyeSamples;

uniform vec3 CameraU;
uniform vec3 CameraV;
uniform vec3 CameraW;
uniform vec3 CameraParams;
uniform vec2 AAOffset;
uniform vec3 CameraPosition;

uniform sampler2D RayDirectionTexture;
uniform sampler2D RayOriginTexture;
uniform sampler2D TexturePolygons;
uniform vec2 PolygonDataStride;

uniform sampler2D TextureBVH;
uniform samplerCube CubeTextureBBoxRootIndices;

uniform sampler2D TextureMaterials;

uniform float MaterialStride;
uniform float MaterialNumRcp;

uniform vec4 OffsetToBBoxMinMax;


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


void next(inout vec2 TriangleIndex, const float offset)
{
	TriangleIndex.y = TriangleIndex.y + PolygonDataStride.y * floor(TriangleIndex.x + PolygonDataStride.x * offset);
	TriangleIndex.x = fract(TriangleIndex.x + PolygonDataStride.x * offset);
}


vec4 LastIntersection;
Intersection raytrace(const Ray ray, const bool cullBackface)
{
	Intersection Result;

	Result.t = 1.0e+30;
	Result.nrm = vec3(0.0);
	Result.pos = vec3(0.0);
	Result.col = vec3(0.0);

	vec3 RayDirection = ray.dir;
	vec3 RayOrigin = ray.org;

	vec4 isect = vec4(-1.0, -1.0, -1.0, 1.0e20);
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

		// read triangle vertices
		vec2 TriangleStartIndex = vec2(BBoxMinTriangleX.w, BBoxMaxTriangleY.w);
		vec2 TriangleIndex = TriangleStartIndex;

		// (px, py, pz, tx)
		// (nx, ny, sgn(nz) * (matid + 1), ty)
		vec3 V0 = texture2D(TexturePolygons, TriangleIndex).xyz; next(TriangleIndex, 1.0); 
		vec3 V1 = texture2D(TexturePolygons, TriangleIndex).xyz; next(TriangleIndex, 1.0); 
		vec3 V2 = texture2D(TexturePolygons, TriangleIndex).xyz; next(TriangleIndex, 1.0); 
		float MaterialIndex = abs(texture2D(TexturePolygons, TriangleIndex).z) - 1.0;

		// perform ray-triangle intersection if it is a leaf node
		if ((BBoxMinTriangleX.w >= 0.0) && BBHit)
		{
			// ray triangle intersection
			//2 cross, 4 dot, 1 /
			//2 cross(12*), 4 dot(12*), 1 /   =  24* 1/
			vec3 p0 = V0;
			vec3 e0 = V1 - V0;
			vec3 e1 = V2 - V0;
			vec3 pv = cross(RayDirection, e1);

			float det = dot(e0, pv);
			if ((cullBackface && (det > 1e-10)) || !cullBackface)
			{
				vec3 tv = RayOrigin - p0;
				vec3 qv = cross(tv, e0);

				vec4 uvt;
				uvt.x = dot(tv, pv);
				uvt.y = dot(RayDirection, qv);
				uvt.z = dot(e1, qv);
				uvt.xyz = uvt.xyz / det;
				uvt.w = 1.0 - uvt.x - uvt.y;

				if (all(greaterThanEqual(uvt, vec4(0.0))) && (uvt.z < isect.a)) 
				{
					Barycentric = uvt.ywx;
					isect = vec4(TriangleStartIndex, MaterialIndex, uvt.z);
				}
			}
			
			/*
			
    // Calculate dimension where the ray direction is maximal.
    ray_coeff_.kz = 0;
    T absDir = std::fabs(ray.dir[0]);
    if (absDir < std::fabs(ray.dir[1])) {
      ray_coeff_.kz = 1;
      absDir = std::fabs(ray.dir[1]);
    }
    if (absDir < std::fabs(ray.dir[2])) {
      ray_coeff_.kz = 2;
      absDir = std::fabs(ray.dir[2]);
    }

    ray_coeff_.kx = ray_coeff_.kz + 1;
    if (ray_coeff_.kx == 3) ray_coeff_.kx = 0;
    ray_coeff_.ky = ray_coeff_.kx + 1;
    if (ray_coeff_.ky == 3) ray_coeff_.ky = 0;

    // Swap kx and ky dimention to preserve widing direction of triangles.
    if (ray.dir[ray_coeff_.kz] < 0.0f) std::swap(ray_coeff_.kx, ray_coeff_.ky);

    // Claculate shear constants.
    ray_coeff_.Sz = 1.0f / ray.dir[ray_coeff_.kz];
    ray_coeff_.Sx = ray.dir[ray_coeff_.kx] * ray_coeff_.Sz;
    ray_coeff_.Sy = ray.dir[ray_coeff_.ky] * ray_coeff_.Sz;
	
	
	const unsigned int f0 = faces_[3 * prim_index + 0];
    const unsigned int f1 = faces_[3 * prim_index + 1];
    const unsigned int f2 = faces_[3 * prim_index + 2];

    const real3<T> p0(get_vertex_addr(vertices_, f0 + 0, vertex_stride_bytes_));
    const real3<T> p1(get_vertex_addr(vertices_, f1 + 0, vertex_stride_bytes_));
    const real3<T> p2(get_vertex_addr(vertices_, f2 + 0, vertex_stride_bytes_));

    const real3<T> A = p0 - ray_org_;
    const real3<T> B = p1 - ray_org_;
    const real3<T> C = p2 - ray_org_;

    const T Ax = A[ray_coeff_.kx] - ray_coeff_.Sx * A[ray_coeff_.kz];
    const T Ay = A[ray_coeff_.ky] - ray_coeff_.Sy * A[ray_coeff_.kz];
    const T Bx = B[ray_coeff_.kx] - ray_coeff_.Sx * B[ray_coeff_.kz];
    const T By = B[ray_coeff_.ky] - ray_coeff_.Sy * B[ray_coeff_.kz];
    const T Cx = C[ray_coeff_.kx] - ray_coeff_.Sx * C[ray_coeff_.kz];
    const T Cy = C[ray_coeff_.ky] - ray_coeff_.Sy * C[ray_coeff_.kz];
	//6*

    T U = Cx * By - Cy * Bx;
    T V = Ax * Cy - Ay * Cx;
    T W = Bx * Ay - By * Ax;
	//6*
	
	//12*

    if (trace_options_.cull_back_face) {
      if (U < 0.0 || V < 0.0 || W < 0.0) return false;
    } else {
      if ((U < 0.0 || V < 0.0 || W < 0.0) && (U > 0.0 || V > 0.0 || W > 0.0)) {
        return false;
      }
    }

    T det = U + V + W;
    if (det == 0.0) return false;

    const T Az = ray_coeff_.Sz * A[ray_coeff_.kz];
    const T Bz = ray_coeff_.Sz * B[ray_coeff_.kz];
    const T Cz = ray_coeff_.Sz * C[ray_coeff_.kz];
    const T D = U * Az + V * Bz + W * Cz;
	//6*
	//18*

    const T rcpDet = 1.0 / det;
    T tt = D * rcpDet;

    if (tt > (*t_inout)) {
      return false;
    }

    (*t_inout) = tt;
    // Use Thomas-Mueller style barycentric coord.
    // U + V + W = 1.0 and interp(p) = U * p0 + V * p1 + W * p2
    // We want interp(p) = (1 - u - v) * p0 + u * v1 + v * p2;
    // => u = V, v = W.
    intersection.u = V * rcpDet;
    intersection.v = W * rcpDet;
	//20*

    return true;
*/
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
	}

	if (isect.x > 0.0)
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


vec3 onb(const vec3 x, const vec3 n)
{
	vec3 u, w, v;
	v = n;
	 
	if (n.z < -0.9999999)
	{
		//u = vec3(0.0, -1.0, 0.0); 
		//w = vec3(-1.0, 0.0, 0.0);
		//return (x.x * u + x.y * v + x.z * w);
		//x.x*(0,-1,0) + x.y*v + x.z*(-1,0,0)
		//(0,-x.x,0) + x.y*v + (-x.z,0,0)
		//x.y*v + (0,-x.x,0) + (-x.z,0,0)
		//x.y*v - (x.z,x.x,0)
		
		//return x.y*v - vec3(x.z,x.x,0); //-2x3*  -6*
		
		u = vec3(x.z, x.x, 0.0);
		
		return x.y*v - u; //-2x3*  -6*
	}
	else
	{
		float a = 1.0 / (1.0 + n.z); 
		float b = -n.x * n.y * a;
		u = vec3(1.0 - n.x * n.x * a, b, -n.x); 
		w = vec3(b, 1.0 - n.y * n.y * a, -n.y);
		return (x.x * u + x.y * v + x.z * w);
	}
}

	// Based on the post http://gpgpu.org/forums/viewtopic.php?t=2591&sid=17051481b9f78fb49fba5b98a5e0f1f3
	// (The page no longer exists as of March 17th, 2015. Please let me know if you see why this code works.)
	const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);
	const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);
	const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);
	const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);

	const vec4 m_2 = vec4(0.5) * m;
	const vec4 m_ = vec4(1) / m;
	
	const vec4 q_ = vec4(1) / q;


float GPURnd(inout vec4 n)
{	
	vec4 beta = floor(n * q_);
	vec4 p = a * (n - beta * q) - beta * r;
	
	beta = (sign(-p)) * m_2 + m_2;
	
	n = (p + beta);

	return fract(dot(n * m_, vec4(1.0, -1.0, 1.0, -1.0)));
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


vec3 glossy_reflect(const vec3 d, const vec3 n, const float g, inout vec4 rndv)
{
	vec3 r = normalize((1.0 - g - g) * reflect(d, n) + (g + g) * n);

	float rnd1 = GPURnd(rndv);
	float rnd2 = GPURnd(rndv);

	//float temp1 = 2.0 * 3.141592 * rnd1;
	float temp1 = 6.283184 * rnd1;
	//float temp2 = sqrt(1.0 - pow(rnd2, 2.0 / (g + 1.0)));
	
	float temp2_ = pow(rnd2, g);
	float temp2 = sqrt(1.0 - temp2_*temp2_);
	vec3 v = vec3(sin(temp1) * temp2, temp2_, cos(temp1) * temp2);

	vec3 result = normalize(onb(v, r));
	if (dot(result, n) < 0.0)
	{
		result = reflect(result, n);
	}

	return result;
}



void main()
{
	vec2 PixelIndex = gl_TexCoord[0].st; 

	vec4 rnd = texture2D(RandomTexture, PixelIndex);
	vec4 rndv = rnd;

	// generate eye ray
	Ray r;
	float FilterResponse;
	{
		// look-at camera
		vec2 PixelPosition = gl_FragCoord.xy + AAOffset - CameraParams.xy;
		vec3 RelativeTargetPosition = PixelPosition.x * CameraU + PixelPosition.y * CameraV + CameraParams.z * CameraW;
		r.dir = normalize(RelativeTargetPosition);
		r.org = CameraPosition;
		FilterResponse = (1.0 - abs(AAOffset.x) / 1.5) * (1.0 - abs(AAOffset.y) / 1.5);

		// thin-lens
		vec3 fp = CameraPosition + r.dir * FocalLength;
		float radius = sqrt(GPURnd(rndv)) * ApertureSize;
		float theta = 6.283184 * GPURnd(rndv);
		vec3 lens = CameraPosition + radius * (CameraU * cos(theta) + CameraV * sin(theta));
		r.dir = normalize(fp - lens);
		r.org = lens;
	}

	vec3 col = vec3(1.0, 1.0, 1.0) * FilterResponse;
	vec3 nrm = vec3(0.0, 0.0, 0.0);
	vec3 pos = vec3(0.0, 0.0, 0.0);
	vec3 emi = vec3(0.0, 0.0, 0.0);
	const float eps = 1e-5;

	for (int j = 0; j < MaxPathLength; j++)
	{
		Intersection i;
		i = raytrace(r, j == 0);
		if (i.t == 1.0e+30)
		{
			nrm = vec3(0.0);
			if (LightSummedArea == 0.0) emi = vec3(0.7, 0.7, 1.0) * col;
			break;
		}

		pos = i.pos;
		nrm = i.nrm;

		if (abs(i.tex.x) < 1e+10) 
		{
			i.col *= texture3D(VolumeTextureTextures, vec3(i.tex, MaterialNumRcp * (i.matid + 0.5))).rgb;
		}

		if (i.brdf == -1)
		{
			nrm = vec3(0.0);
			emi = col * i.col;
			break;
		}
		else if (i.brdf == 0)
		{
			col = col * i.col;
			break;
		}
		else if ((i.brdf == 1) || (i.brdf == 4))
		{
			r.org = pos + eps * i.gnrm;
			
			//if (i.brdf == 4) nrm = glossy_reflect(nrm, nrm, 1.0 / pow((1.0 - i.g) * 0.5, 2.71828), rndv);
			//tigra: 1/ 1* => 1/
			//if (i.brdf == 4) nrm = glossy_reflect(nrm, nrm, pow(2.0 / (1.0 - i.g), 2.71828), rndv);
			if (i.brdf == 4) nrm = glossy_reflect(nrm, nrm, i.g, rndv);
			r.dir = reflect(r.dir, nrm);

			if (dot(r.dir, i.gnrm) < 0.0) r.dir = -r.dir;
			col = col * i.col;
		}
		else if ((i.brdf == 2) || (i.brdf == 5))
		{
			//if (i.brdf == 5) nrm = glossy_reflect(nrm, nrm, 1.0 / pow((1.0 - i.g) * 0.5, 2.71828), rndv);
			//if (i.brdf == 5) nrm = glossy_reflect(nrm, nrm, pow(2.0 / (1.0 - i.g), 2.71828), rndv);
			if (i.brdf == 5) nrm = glossy_reflect(nrm, nrm, i.g, rndv);

			// specular refraction
			float ln = dot(nrm, r.dir);
			float eta = i.eta;
			if (ln < 0.0)
			{
				// in
				float Re = Fresnel(-r.dir, nrm, 1.0, eta);
				if (GPURnd(rndv) < Re)
				{
					r.org = pos + eps * i.gnrm;
					r.dir = reflect(r.dir, nrm);
					if (dot(r.dir, i.gnrm) < 0.0) r.dir = -r.dir;
				}
				else
				{
					col = col * i.col;
					r.org = pos - eps * i.gnrm;
					r.dir = refract(r.dir, nrm, 1.0 / eta);
					if (dot(r.dir, i.gnrm) > 0.0) r.dir = -r.dir;
				}
			}
			else
			{
				// out
				float Re = Fresnel(-r.dir, -nrm, eta, 1.0);
				col = col * i.col;
				if (GPURnd(rndv) < Re)
				{
					r.org = pos - eps * i.gnrm;
					r.dir = reflect(r.dir, -nrm);
					if (dot(r.dir, -i.gnrm) < 0.0) r.dir = -r.dir;
				}
				else
				{
					r.org = pos + eps * i.gnrm;
					r.dir = refract(r.dir, -nrm, eta);
					if (dot(r.dir, -i.gnrm) > 0.0) r.dir = -r.dir;
				}
			}
		}
		else if ((i.brdf == 3) || (i.brdf == 6))
		{
			//if (i.brdf == 6) nrm = glossy_reflect(nrm, nrm, 1.0 / pow((1.0 - i.g) * 0.5, 2.71828), rndv);
			//if (i.brdf == 6) nrm = glossy_reflect(nrm, nrm, pow(2.0 / (1.0 - i.g), 2.71828), rndv);
			if (i.brdf == 6) nrm = glossy_reflect(nrm, nrm, i.g, rndv);

			// specular reflection
			float ln = -abs(dot(nrm, r.dir));
			float eta = i.eta;
			float Re = Fresnel(-r.dir, nrm, 1.0, eta);
			if (GPURnd(rndv) < Re)
			{
				r.org = pos + eps * i.gnrm;
				r.dir = reflect(r.dir, nrm);
				if (dot(r.dir, i.gnrm) < 0.0) r.dir = -r.dir;
			}
			else
			{
				col = col * i.col;
				break;
			}
		}
	}

	// emission
	vec4 QueryEmissionPhotonCount = texture2D(QueryEmissionPhotonCountTexture, PixelIndex);
	float a = 1.0 / float(NumEyeSamples);
	QueryEmissionPhotonCount.rgb = QueryEmissionPhotonCount.rgb * (1.0 - a) + emi * a;
	 
	// eye ray tracing
	gl_FragData[0] = vec4(pos, r.dir.r);
	gl_FragData[1] = vec4(col, r.dir.g);
	gl_FragData[2] = vec4(nrm, r.dir.b);
	gl_FragData[3] = rndv;
	gl_FragData[4] = LastIntersection;
	gl_FragData[5] = QueryEmissionPhotonCount;
}
