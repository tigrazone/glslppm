#ifndef _GPURT_H
#define _GPURT_H

#include <stdlib.h>

#ifdef __APPLE__
	#include <GLUT/glut.h>
#else
	#include "GLee.h"
	#include "gl/glut.h"
#endif

#include "vector.h"
#include "mesh.h"
#include "bvh.h"
#include "camera.h"


// temporary data structure for creating textures for polygons
struct TBBoxGPU
{
	TVector3 min, max;
	TVector2 hit, miss;
	int tri;
};


// GPURT class
class CGPURT
{
	public:
		CGPURT()
		{
			textureTriangles = 0;
			textureBVHs = 0;

			TextureMaterials = 0;
			TextureLightSources = 0;
			VolumeTextureTextures = 0;

			CubeTextureBBoxRootIndices = 0;
		}

		// textures
		GLuint textureTriangles;
		GLuint textureBVHs;
		GLuint VolumeTextureTextures;

		GLuint CubeTextureBBoxRootIndices;
		GLuint TextureMaterials;
		float MaterialDataStride;
		GLuint TextureLightSources;
		TVector2 PolygonDataStride;
		int BBoxDataSizeX, BBoxDataSizeY;

		TMesh mesh;
		TCamera camera;

		void Release();
		void PrecalculateMeshData();

	private:
		void BuildBVH();
		void CreateTextures();
		GLuint CreateTexture(const int TextureIndex, GLenum format, const int Width, const int Height);
		TBVH bvh;
};

#endif