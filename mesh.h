#ifndef _MESH_H
#define _MESH_H


#include "vector.h"
#include <vector>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bbox.h"


class TTriangle
{
	public:
		TVector3 positions[3];
		TVector3 normals[3];
		TVector2 texcoords[3];

		TBBox bbox;
		TVector3 centroid;
		int idMaterial;

	private:
};


class TMaterial
{
	public:
		std::string name;
		TVector3 color;
		int brdf;
		float eta;
		float specularity;

		TVector3 Ka, Kd, Ks;
		float Ns;

		bool isTextured;
		std::vector<unsigned char> texture;
		unsigned int textureWidth;
		unsigned int textureHeight;
};


class TMesh
{
	public:
		std::vector<TTriangle> triangles;
		std::vector<TMaterial> materials;

		std::vector<float> lightsCDF;
		std::vector<int> lightsIndices;
		float lightsArea;

		TBBox bbox;

		void LoadOBJ(char* FileName, TVector3 Position, float Scale);
		void CalculateBBox();
		void PrepareLightSources();
		void Release();

	private:
};


#endif
