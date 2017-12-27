// glslppm, a GPU implementation of Stochastic Progressive Photon Mapping by T. Hachisuka, 2015
// March 17, 2015: Initial release.

#include <stdlib.h>
#include "gpurt.h"

#ifdef __APPLE__
	#include <GLUT/glut.h>
#else
	#include "GLee.h"
	#include "gl/glut.h"
#endif

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>


static TVector3 CanonicalCameraPosition;
static float FieldOfView;
static TVector3 LookAtPosition;

static CGPURT gpurt;

GLuint PSDraw, PSHash, PSEyeRayTrace, PSPhotonTrace, PSProgressiveUpdate, PSRadianceEstimate, PVSScatter, PVSCorrection, PSMax, PSMin, PSSum;

static float NumPhotons = 0.0f;
static unsigned int NumEyeSamples = 0;
static unsigned int FrameCount = 0;
static float FocalLength = 13.0f;
static float ApertureSize = 0.0f;
static int m_nextUpdate = 0;
static int startedTime = 0;

const unsigned int MaxNumberOfBounces = 10;
const unsigned int ImageResolution = 512;
const unsigned int PhotonBufferSize = ImageResolution;
const unsigned int HashResolution = ImageResolution;
const float InitialFootprint = 2.5f;

GLuint QueryPositionTexture;
GLuint QueryNormalTexture;
GLuint QueryEmissionPhotonCountTexture;
GLuint QueryFluxRadiusTexture;
GLuint QueryReflectanceTexture;
GLuint QueryIntersectionTexture;

GLuint PhotonIndexTexture;
GLuint PhotonFluxTexture;
GLuint PhotonPositionTexture;
GLuint PhotonDirectionTexture;
GLuint PhotonHashTexture;
GLuint PhotonCorrectionTexture;
GLuint RandomPhotonTexture;
GLuint RandomEyeRayTexture;
GLuint PhotonIntersectionTexture;
GLuint PhotonEmittedFlagTexture;


GLuint EyeRayTraceSurface;
GLuint PhotonRayTraceSurface;
GLuint PhotonIndexSurface;
GLuint QueryPointSurface;
GLuint PhotonHashSurface;
GLuint PhotonHashDepthBuffer;
GLuint PhotonCorrectionSurface;

GLuint MinMaxAveTextureQuery;
GLuint MinMaxAveSurfaceQuery;
GLuint MinMaxAveTexturePhoton;
GLuint MinMaxAveSurfacePhoton;

GLuint FragmentsVBO;


namespace XORShift
{
	// XOR shift PRNG
	static unsigned int m_x = 123456789;
	static unsigned int m_y = 362436069;
	static unsigned int m_z = 521288629;
	static unsigned int m_w = 88675123;

	inline float m_frand()
	{
		const unsigned int t = m_x ^ (m_x << 11);
		m_x = m_y; m_y = m_z; m_z = m_w;
		return (m_w = (m_w ^ (m_w >> 19)) ^ (t ^ (t >> 8))) * (1.0f / 4294967295.0f);
	}
}


void m_SetTexture(const int TextureIndex, GLuint Texture)
{
	glActiveTexture(GL_TEXTURE0 + TextureIndex);
	glBindTexture(GL_TEXTURE_2D, Texture);
}

void m_SetCubeTexture(const int CubeTextureIndex, GLuint CubeTexture)
{
	glActiveTexture(GL_TEXTURE0 + CubeTextureIndex);
	glBindTexture(GL_TEXTURE_CUBE_MAP, CubeTexture);
}

void m_SetVolumeTexture(const int VolumeTextureIndex, GLuint VolumeTexture)
{
	glActiveTexture(GL_TEXTURE0 + VolumeTextureIndex);
	glBindTexture(GL_TEXTURE_3D, VolumeTexture);
}


void RandomizeTextures()
{
	std::vector<TVector4> TempData;

	// for eye rays
	TempData.resize(ImageResolution * ImageResolution);
	m_SetTexture(0, RandomEyeRayTexture);
	for (int j = 0; j < ImageResolution; j++)
	{
		for (int i = 0; i < ImageResolution; i++)
		{
			TempData[i + j * ImageResolution] = TVector4(XORShift::m_frand() * 4194304.0, XORShift::m_frand() * 4194304.0, XORShift::m_frand() * 4194304.0, XORShift::m_frand() * 4194304.0);
		}
	}
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, ImageResolution, ImageResolution, GL_RGBA, GL_FLOAT, &TempData[0]);

	// for photons
	TempData.resize(PhotonBufferSize * PhotonBufferSize);
	m_SetTexture(0, RandomPhotonTexture);
	for (int j = 0; j < PhotonBufferSize; j++)
	{
		for (int i = 0; i < PhotonBufferSize; i++)
		{
			TempData[i + j * PhotonBufferSize] = TVector4(XORShift::m_frand() * 4194304.0, XORShift::m_frand() * 4194304.0, XORShift::m_frand() * 4194304.0, XORShift::m_frand() * 4194304.0);
		}
	}
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, PhotonBufferSize, PhotonBufferSize, GL_RGBA, GL_FLOAT, &TempData[0]);
}


void drawQuad(const int w, const int h)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, w, 0.0, h);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glViewport(0, 0, w, h);

	glBegin(GL_QUADS);
		glTexCoord2f(0.0, 0.0); glVertex2f(0.0, 0.0);
		glTexCoord2f(1.0, 0.0); glVertex2f(  w, 0.0);
		glTexCoord2f(1.0, 1.0); glVertex2f(  w,   h);
		glTexCoord2f(0.0, 1.0); glVertex2f(0.0,   h);
	glEnd();
}

void drawQuadwithTex(const int w, const int h, const float s, const float t)
{
	glBegin(GL_QUADS);
		glTexCoord2f(0.0, 0.0); glVertex2f(0.0, 0.0);
		glTexCoord2f(  s, 0.0); glVertex2f(  w, 0.0);
		glTexCoord2f(  s,   t); glVertex2f(  w,   h);
		glTexCoord2f(0.0,   t); glVertex2f(0.0,   h);
	glEnd();
}

TVector4 ReduceTexture(const GLuint MinMaxAveSurface, const GLuint MinMaxAveTexture, const GLuint Texture, const GLuint Shader, const unsigned int Resolution)
{
	// this function assumes ImageResolution = 2^t and the image is a square
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, MinMaxAveSurface);
	glUseProgram(Shader);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, Resolution, 0.0, Resolution);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glViewport(0, 0, Resolution, Resolution);

	// first pass reduces and copies Texture into MinMaxAveTexture
	unsigned int Level = 1;
	unsigned int ReducedBufferSize = Resolution >> Level;
	float TextureOffset = 1.0f / float(1 << Level);

	glUniform1i(glGetUniformLocation(Shader, "Texture"), 15); m_SetTexture(15, Texture);
	glUniform2f(glGetUniformLocation(Shader, "Offset"), TextureOffset, TextureOffset);
	drawQuadwithTex(ReducedBufferSize, ReducedBufferSize, TextureOffset, TextureOffset);

	// remaining passes keep reducing MinMaxAveTexture
	unsigned int numPasses = (unsigned)(log((double)(Resolution >> Level)) / log(2.0)) + 1;
	TVector4 result;
	for (unsigned int i = 0; i < numPasses; i++)
	{
		Level++;
		TextureOffset = 1.0f / float(1 << Level);
		ReducedBufferSize = Resolution >> Level;

		glUniform1i(glGetUniformLocation(Shader, "Texture"), 15); m_SetTexture(15, MinMaxAveTexture);
		glUniform2f(glGetUniformLocation(Shader, "Offset"), TextureOffset, TextureOffset);
		drawQuadwithTex(ReducedBufferSize, ReducedBufferSize, TextureOffset, TextureOffset);

		// make sure that the rendering process is done
		glReadPixels(0, 0, 1, 1, GL_RGBA, GL_FLOAT, &result);
	}
	return result;
}


static TVector4 BBMax, BBMin;
static float InitialRadius, GridScale;
void m_display(void)
{
	TVector4 BBoxOffsets;
	BBoxOffsets.x = (float)(                             0.5f) / (float)(gpurt.BBoxDataSizeX * 2.0f);
	BBoxOffsets.y = (float)(gpurt.BBoxDataSizeY * 3.0f + 0.5f) / (float)(gpurt.BBoxDataSizeY * 4.0f);
	BBoxOffsets.z = (float)(gpurt.BBoxDataSizeX        + 0.5f) / (float)(gpurt.BBoxDataSizeX * 2.0f);
	BBoxOffsets.w = (float)(gpurt.BBoxDataSizeY * 3.0f + 0.5f) / (float)(gpurt.BBoxDataSizeY * 4.0f);

	if (FrameCount == 0)
	{
		// emission & local photon count
		std::vector<TVector4> TempData(ImageResolution * ImageResolution);
		m_SetTexture(5, QueryEmissionPhotonCountTexture);
		for (int j = 0; j < ImageResolution; j++)
		{
			for (int i = 0; i < ImageResolution; i++)
			{
				TempData[i + j * ImageResolution] = TVector4(0.0f, 0.0f, 0.0f, 0.0f);
			}
		}
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, ImageResolution, ImageResolution, GL_RGBA, GL_FLOAT, &TempData[0]);
	}

	// balance the cost of eye ray tracing and photon ray tracing
	if ((FrameCount % MaxNumberOfBounces) == 0)
	{
		// eye ray tracing
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, EyeRayTraceSurface);

		glUseProgram(PSEyeRayTrace);

		// ray tracing parameters
		glUniform4f(glGetUniformLocation(PSEyeRayTrace, "OffsetToBBoxMinMax"), BBoxOffsets.x, BBoxOffsets.y, BBoxOffsets.z, BBoxOffsets.w);
		glUniform1i(glGetUniformLocation(PSEyeRayTrace, "TexturePolygons"), 2); m_SetTexture(2, gpurt.textureTriangles);
		glUniform1i(glGetUniformLocation(PSEyeRayTrace, "CubeTextureBBoxRootIndices"), 5); m_SetCubeTexture(5, gpurt.CubeTextureBBoxRootIndices);
		glUniform1i(glGetUniformLocation(PSEyeRayTrace, "TextureBVH"), 6); m_SetTexture(6, gpurt.textureBVHs);

		glUniform1i(glGetUniformLocation(PSEyeRayTrace, "VolumeTextureTextures"), 11); m_SetVolumeTexture(11, gpurt.VolumeTextureTextures);

		// material data
		glUniform1i(glGetUniformLocation(PSEyeRayTrace, "TextureMaterials"), 10); m_SetTexture(10, gpurt.TextureMaterials);
		glUniform1f(glGetUniformLocation(PSEyeRayTrace, "MaterialStride"), gpurt.MaterialDataStride);
		glUniform1f(glGetUniformLocation(PSEyeRayTrace, "MaterialNumRcp"), 1.0f / float(gpurt.mesh.materials.size()));
		glUniform1f(glGetUniformLocation(PSEyeRayTrace, "LightSummedArea"), gpurt.mesh.lightsArea);

		// camera parameters
		TVector4 VData;
		VData.x = gpurt.camera.origin.x;
		VData.y = gpurt.camera.origin.y;
		VData.z = gpurt.camera.origin.z;
		glUniform3f(glGetUniformLocation(PSEyeRayTrace, "CameraPosition"), VData.x, VData.y, VData.z);

		glUniform3f(glGetUniformLocation(PSEyeRayTrace, "CameraU"), gpurt.camera.u.x, gpurt.camera.u.y, gpurt.camera.u.z);
		glUniform3f(glGetUniformLocation(PSEyeRayTrace, "CameraV"), gpurt.camera.v.x, gpurt.camera.v.y, gpurt.camera.v.z);
		glUniform3f(glGetUniformLocation(PSEyeRayTrace, "CameraW"), gpurt.camera.w.x, gpurt.camera.w.y, gpurt.camera.w.z);

		VData.x = (float)gpurt.camera.width * 0.5f;
		VData.y = (float)gpurt.camera.height * 0.5f;
		VData.z = gpurt.camera.distance;
		glUniform3f(glGetUniformLocation(PSEyeRayTrace, "CameraParams"), VData.x, VData.y, VData.z);

		// antialiasing offset
		VData.x = (XORShift::m_frand() - 0.5f) * 1.25f;
		VData.y = (XORShift::m_frand() - 0.5f) * 1.25f;
		glUniform2f(glGetUniformLocation(PSEyeRayTrace, "AAOffset"), VData.x, VData.y);

		// some extra parameters
		glUniform1i(glGetUniformLocation(PSEyeRayTrace, "RandomTexture"), 0); m_SetTexture(0, RandomEyeRayTexture);
		glUniform1i(glGetUniformLocation(PSEyeRayTrace, "QueryEmissionPhotonCountTexture"), 1); m_SetTexture(1, QueryEmissionPhotonCountTexture);
		glUniform1f(glGetUniformLocation(PSEyeRayTrace, "FocalLength"), FocalLength);
		glUniform1i(glGetUniformLocation(PSEyeRayTrace, "MaxPathLength"), MaxNumberOfBounces);
		glUniform1f(glGetUniformLocation(PSEyeRayTrace, "ApertureSize"), ApertureSize);
		glUniform2f(glGetUniformLocation(PSEyeRayTrace, "PolygonDataStride"), gpurt.PolygonDataStride.x, gpurt.PolygonDataStride.y);
		glUniform1i(glGetUniformLocation(PSEyeRayTrace, "NumEyeSamples"), NumEyeSamples + 1);

		NumEyeSamples++;
		drawQuad(ImageResolution, ImageResolution);

		BBMax = ReduceTexture(MinMaxAveSurfaceQuery, MinMaxAveTextureQuery, QueryPositionTexture, PSMax, ImageResolution);
		BBMin = ReduceTexture(MinMaxAveSurfaceQuery, MinMaxAveTextureQuery, QueryPositionTexture, PSMin, ImageResolution);
		float BBSize = 0.0f;
		for (int i = 0; i < 3; i++)
		{
			BBSize += BBMax[i] - BBMin[i];
		}

		// initial radius estimation
		InitialRadius = (BBSize / 3.0f) / (float)(ImageResolution) * InitialFootprint;

		// expand the bounding box
		for (int i = 0; i < 3; i++)
		{
			BBMin[i] -= InitialRadius;
			BBMax[i] += InitialRadius;
		}

		// hashed grid resolution
		GridScale = 0.5f / InitialRadius;
	}


	// initialized photons
	if (FrameCount == 0)
	{
		std::vector<TVector4> TempData(ImageResolution * ImageResolution);
		for (int j = 0; j < ImageResolution; j++)
		{
			for (int i = 0; i < ImageResolution; i++)
			{
				TempData[i + j * ImageResolution] = TVector4(0.0f, 0.0f, 0.0f, 0.0f);
			}
		}

		// accumulated (unnormalized) flux & radius
		m_SetTexture(6, QueryFluxRadiusTexture);
		for (int j = 0; j < ImageResolution; j++)
		{
			for (int i = 0; i < ImageResolution; i++)
			{
				TempData[i + j * ImageResolution] = TVector4(0.0f, 0.0f, 0.0f, InitialRadius);
			}
		}
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, ImageResolution, ImageResolution, GL_RGBA, GL_FLOAT, &TempData[0]);

		// photon intersection
		TempData.resize(PhotonBufferSize * PhotonBufferSize);
		m_SetTexture(7, PhotonIntersectionTexture);
		for (int j = 0; j < PhotonBufferSize; j++)
		{
			for (int i = 0; i < PhotonBufferSize; i++)
			{
				TempData[i + j * PhotonBufferSize] = TVector4(-1.0f, -1.0f, 0.0f, 1.0e+30f);
			}
		}
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, PhotonBufferSize, PhotonBufferSize, GL_RGBA, GL_FLOAT, &TempData[0]);

		TempData.clear();
		NumPhotons = 0.0f;
	}



	{
		// photon tracing
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, PhotonRayTraceSurface);
			glUseProgram(PSPhotonTrace);

			// ray tracing
			glUniform4f(glGetUniformLocation(PSPhotonTrace, "OffsetToBBoxMinMax"), BBoxOffsets.x, BBoxOffsets.y, BBoxOffsets.z, BBoxOffsets.w);
			glUniform1i(glGetUniformLocation(PSPhotonTrace, "TexturePolygons"), 2); m_SetTexture(2, gpurt.textureTriangles);
			glUniform1i(glGetUniformLocation(PSPhotonTrace, "CubeTextureBBoxRootIndices"), 5); m_SetCubeTexture(5, gpurt.CubeTextureBBoxRootIndices);
			glUniform1i(glGetUniformLocation(PSPhotonTrace, "TextureBVH"), 6); m_SetTexture(6, gpurt.textureBVHs);
			glUniform2f(glGetUniformLocation(PSPhotonTrace, "PolygonDataStride"), gpurt.PolygonDataStride.x, gpurt.PolygonDataStride.y);

			glUniform1i(glGetUniformLocation(PSPhotonTrace, "PhotonIntersectionTexture"), 12); m_SetTexture(12, PhotonIntersectionTexture);
			glUniform1i(glGetUniformLocation(PSPhotonTrace, "PhotonPositionTexture"), 13); m_SetTexture(13, PhotonPositionTexture);
			glUniform1i(glGetUniformLocation(PSPhotonTrace, "PhotonFluxTexture"), 14); m_SetTexture(14, PhotonFluxTexture);
			glUniform1i(glGetUniformLocation(PSPhotonTrace, "PhotonDirectionTexture"), 15); m_SetTexture(15, PhotonDirectionTexture);

			// brdfs
			glUniform1i(glGetUniformLocation(PSPhotonTrace, "TextureMaterials"), 10); m_SetTexture(10, gpurt.TextureMaterials);
			glUniform1f(glGetUniformLocation(PSPhotonTrace, "MaterialStride"), gpurt.MaterialDataStride);
			glUniform1f(glGetUniformLocation(PSPhotonTrace, "MaterialNumRcp"), 1.0f / float(gpurt.mesh.materials.size()));

			// material data
			glUniform1i(glGetUniformLocation(PSPhotonTrace, "VolumeTextureTextures"), 9); m_SetVolumeTexture(9, gpurt.VolumeTextureTextures);
			glUniform1i(glGetUniformLocation(PSPhotonTrace, "TextureLightSources"), 11); m_SetTexture(11, gpurt.TextureLightSources);
			glUniform1f(glGetUniformLocation(PSPhotonTrace, "LightSourceStride"), 1.0f / (float)gpurt.mesh.lightsCDF.size());
			glUniform1f(glGetUniformLocation(PSPhotonTrace, "LightSummedArea"), gpurt.mesh.lightsArea);

			glUniform1i(glGetUniformLocation(PSPhotonTrace, "RandomTexture"), 0); m_SetTexture(0, RandomPhotonTexture);
			glUniform1i(glGetUniformLocation(PSPhotonTrace, "MaxPathLength"), MaxNumberOfBounces);

			// for IBL
			TVector4 SceneBSphere;
			SceneBSphere.x = (gpurt.mesh.bbox.max.x + gpurt.mesh.bbox.min.x) * 0.5f;
			SceneBSphere.y = (gpurt.mesh.bbox.max.y + gpurt.mesh.bbox.min.y) * 0.5f;
			SceneBSphere.z = (gpurt.mesh.bbox.max.z + gpurt.mesh.bbox.min.z) * 0.5f;
			SceneBSphere.w = (gpurt.mesh.bbox.max - gpurt.mesh.bbox.min).length() * 0.5f;
			glUniform4f(glGetUniformLocation(PSPhotonTrace, "SceneBSphere"), SceneBSphere.x, SceneBSphere.y, SceneBSphere.z, SceneBSphere.w);

		drawQuad(PhotonBufferSize, PhotonBufferSize);
		TVector4 NumCurrentEmittedPhotons = ReduceTexture(MinMaxAveSurfacePhoton, MinMaxAveTexturePhoton, PhotonEmittedFlagTexture, PSSum, PhotonBufferSize);
		NumPhotons += int(NumCurrentEmittedPhotons.x);
	}


	{
		// build a stochastic hashed grid

		// compute the hash values of the photons
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, PhotonIndexSurface);
			glUseProgram(PSHash);

			glUniform4f(glGetUniformLocation(PSHash, "BufInfo"), HashResolution, HashResolution, 1.0f / float(HashResolution), 1.0f / float(HashResolution));
			glUniform1i(glGetUniformLocation(PSHash, "HashNum"), HashResolution * HashResolution);
			glUniform1f(glGetUniformLocation(PSHash, "GridScale"), GridScale);
			glUniform1f(glGetUniformLocation(PSHash, "HashScale1"), 4194304.0f / GridScale);
			glUniform3f(glGetUniformLocation(PSHash, "BBoxMin"), BBMin.x, BBMin.y, BBMin.z);
			glUniform1i(glGetUniformLocation(PSHash, "PhotonPositionTexture"), 10); m_SetTexture(10, PhotonPositionTexture);
			glUniform1i(glGetUniformLocation(PSHash, "PhotonFluxTexture"), 11); m_SetTexture(11, PhotonFluxTexture);
		drawQuad(PhotonBufferSize, PhotonBufferSize);


		// random write photons into the hashed buffer
		glEnable(GL_DEPTH_TEST);
			glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, PhotonHashSurface);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

				glUseProgram(PVSScatter);

				glUniform4f(glGetUniformLocation(PVSScatter, "BufInfo") , HashResolution, HashResolution, 1.0f / float(HashResolution), 1.0f / float(HashResolution));
				glUniform1i(glGetUniformLocation(PVSScatter, "PhotonIndexTexture"), 8); m_SetTexture(8, PhotonIndexTexture);
				glUniform1f(glGetUniformLocation(PVSScatter, "PhotonBufferSize1"), 1.0f / float(PhotonBufferSize) );

				glMatrixMode(GL_PROJECTION);
					glLoadIdentity();
					gluOrtho2D(0.0, HashResolution, 0.0, HashResolution);
				glMatrixMode(GL_MODELVIEW);
					glLoadIdentity();
				glViewport(0, 0, HashResolution, HashResolution);

			glBindBufferARB(GL_ARRAY_BUFFER_ARB, FragmentsVBO);
				glEnableClientState(GL_VERTEX_ARRAY);
					glVertexPointer(2, GL_FLOAT, 0, 0);
					glDrawArrays(GL_POINTS, 0, PhotonBufferSize * PhotonBufferSize);
				glDisableClientState(GL_VERTEX_ARRAY);
			glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
		glDisable(GL_DEPTH_TEST);


		// count the number of overlapped photons in the hashed grid
		// - this is necessary to make the estimation unbiased (essentially the Russian roulette technique)
		glEnable(GL_BLEND);
			glBlendFunc(GL_ONE, GL_ONE);
			glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, PhotonCorrectionSurface);
				glClear(GL_COLOR_BUFFER_BIT);

				glUseProgram(PVSCorrection);

				glUniform4f(glGetUniformLocation(PVSCorrection, "BufInfo"), HashResolution, HashResolution, 1.0f / float(HashResolution), 1.0f / float(HashResolution));
				glUniform1i(glGetUniformLocation(PVSCorrection, "PhotonIndexTexture"), 8); 	m_SetTexture(8, PhotonIndexTexture);

				glMatrixMode(GL_PROJECTION);
					glLoadIdentity();
					gluOrtho2D(0.0, HashResolution, 0.0, HashResolution);
				glMatrixMode(GL_MODELVIEW);
					glLoadIdentity();
				glViewport(0, 0, HashResolution, HashResolution);

			glBindBufferARB(GL_ARRAY_BUFFER_ARB, FragmentsVBO);
				glEnableClientState(GL_VERTEX_ARRAY);
					glVertexPointer(2, GL_FLOAT, 0, 0);
					glDrawArrays(GL_POINTS, 0, PhotonBufferSize * PhotonBufferSize);
				glDisableClientState(GL_VERTEX_ARRAY);
			glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
		glDisable(GL_BLEND);


		// radiance estimation
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, QueryPointSurface);
			glUseProgram(PSProgressiveUpdate);

			// the maximum hash index
			TVector3 HashMax;
			HashMax.x = fabs(BBMax[0] + InitialRadius - BBMin[0]) * GridScale;
			HashMax.y = fabs(BBMax[1] + InitialRadius - BBMin[1]) * GridScale;
			HashMax.z = fabs(BBMax[2] + InitialRadius - BBMin[2]) * GridScale;
			glUniform3f(glGetUniformLocation(PSProgressiveUpdate, "HashMax"), HashMax.x, HashMax.y, HashMax.z);

			glUniform4f(glGetUniformLocation(PSProgressiveUpdate, "BufInfo"), HashResolution, HashResolution, 1.0f / float(HashResolution), 1.0f / float(HashResolution));
			glUniform1f(glGetUniformLocation(PSProgressiveUpdate, "HashNum"), HashResolution * HashResolution);
			glUniform1f(glGetUniformLocation(PSProgressiveUpdate, "GridScale"), GridScale);
			glUniform1f(glGetUniformLocation(PSProgressiveUpdate, "HashScale1"), 4194304.0f / GridScale);
			glUniform1f(glGetUniformLocation(PSProgressiveUpdate, "Alpha"), 0.7f);
			glUniform3f(glGetUniformLocation(PSProgressiveUpdate, "BBoxMin"), BBMin.x, BBMin.y, BBMin.z);

			glUniform1i(glGetUniformLocation(PSProgressiveUpdate, "HashedPhotonTexture"), 0); m_SetTexture(0, PhotonHashTexture);
			glUniform1i(glGetUniformLocation(PSProgressiveUpdate, "PhotonCorrectionTexture"), 1); m_SetTexture(1, PhotonCorrectionTexture);
			glUniform1i(glGetUniformLocation(PSProgressiveUpdate, "QueryNormalTexture"), 2); m_SetTexture(2, QueryNormalTexture);
			glUniform1i(glGetUniformLocation(PSProgressiveUpdate, "QueryPositionTexture"), 3); m_SetTexture(3, QueryPositionTexture);
			glUniform1i(glGetUniformLocation(PSProgressiveUpdate, "QueryEmissionPhotonCountTexture"), 5); m_SetTexture(5, QueryEmissionPhotonCountTexture);
			glUniform1i(glGetUniformLocation(PSProgressiveUpdate, "QueryFluxRadiusTexture"), 6); m_SetTexture(6, QueryFluxRadiusTexture);
			glUniform1i(glGetUniformLocation(PSProgressiveUpdate, "QueryReflectanceTexture"), 7); m_SetTexture(7, QueryReflectanceTexture);
			glUniform1i(glGetUniformLocation(PSProgressiveUpdate, "PhotonFluxTexture"), 9); m_SetTexture(9, PhotonFluxTexture);
			glUniform1i(glGetUniformLocation(PSProgressiveUpdate, "PhotonPositionTexture"), 10); m_SetTexture(10, PhotonPositionTexture);
			glUniform1i(glGetUniformLocation(PSProgressiveUpdate, "PhotonDirectionTexture"), 11); m_SetTexture(11, PhotonDirectionTexture);
			glUniform1i(glGetUniformLocation(PSProgressiveUpdate, "QueryIntersectionTexture"), 14); m_SetTexture(14, QueryIntersectionTexture);
		drawQuad(ImageResolution, ImageResolution);
	}


	// rendering
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
		glUseProgram(PSRadianceEstimate);

		glUniform1f(glGetUniformLocation(PSRadianceEstimate, "TotalPhotonNum"), NumPhotons);
		glUniform1i(glGetUniformLocation(PSRadianceEstimate, "QueryEmissionPhotonCountTexture"), 5); m_SetTexture(5, QueryEmissionPhotonCountTexture);
		glUniform1i(glGetUniformLocation(PSRadianceEstimate, "QueryFluxRadiusTexture"), 6); m_SetTexture(6, QueryFluxRadiusTexture);
	drawQuad(ImageResolution, ImageResolution);

/*
	// debug output
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
		glClear(GL_COLOR_BUFFER_BIT);
		glUseProgram(PSDraw);
		m_SetTexture(0, PhotonHashTexture);
		glUniform1i(glGetUniformLocation(PSDraw, "input_tex"), 0);
	drawQuad(ImageResolution, ImageResolution);
//*/

	// update
	glFinish();
	FrameCount++;
	if ((glutGet(GLUT_ELAPSED_TIME) - m_nextUpdate) > 0)
	{
		std::stringstream s;
		float NumMPaths = (NumPhotons + NumEyeSamples * ImageResolution * ImageResolution) / (1024.0 * 1024.0);
		s << ((glutGet(GLUT_ELAPSED_TIME) - startedTime) / 1000) << "sec, " << (NumMPaths / float((glutGet(GLUT_ELAPSED_TIME) - startedTime) / 1000)) << "M paths/sec, " << NumMPaths << "M paths" << std::endl;
		glutSetWindowTitle(s.str().c_str());
		m_nextUpdate = glutGet(GLUT_ELAPSED_TIME) + 1000;
	}
}


// error checking for GLSL (from http://www.mathematik.uni-dortmund.de/~goeddeke/gpgpu/tutorial.html)
void m_printInfoLogs(GLuint obj, GLuint shader)
{
	int infologLength = 0;
	int charsWritten  = 0;
	char *infoLog;

	glGetProgramiv(obj, GL_INFO_LOG_LENGTH, &infologLength);
	if (infologLength > 1)
	{
		infoLog = (char *)malloc(infologLength);
		glGetProgramInfoLog(obj, infologLength, &charsWritten, infoLog);
		std::cerr << infoLog << std::endl;
		free(infoLog);
	}
	glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infologLength);
	if (infologLength > 1)
	{
		infoLog = (char *)malloc(infologLength);
		glGetShaderInfoLog(shader, infologLength, &charsWritten, infoLog);
		std::cerr << infoLog << std::endl;
		free(infoLog);
	}
}


GLuint m_CreateFullShader(const char* vertex_shader_path, const char* fragment_shader_path)
{
	// create a fragment shader and a vertex shader
	GLuint p = glCreateProgram();
	std::cerr << "compiling " << vertex_shader_path << "...";
	{
		std::ifstream ifs(vertex_shader_path, std::ios::binary);
		ifs.seekg(0, std::ios::end);
		std::ifstream::pos_type filesize = ifs.tellg();
		ifs.seekg(0, std::ios::beg);
		std::vector<char> bytes(filesize);
		ifs.read(&bytes[0], filesize);
		GLint size = filesize;
		const char* c = &bytes[0];

		GLuint s = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(s,1, &c, &size);
		glCompileShader(s);
		glAttachShader(p,s);
		m_printInfoLogs(p, s);
	}
	std::cerr << "done." << std::endl;
	std::cerr << "compiling " << fragment_shader_path << "...";
	{
		std::ifstream ifs(fragment_shader_path, std::ios::binary);
		ifs.seekg(0, std::ios::end);
		std::ifstream::pos_type filesize = ifs.tellg();
		ifs.seekg(0, std::ios::beg);
		std::vector<char> bytes(filesize);
		ifs.read(&bytes[0], filesize);
		GLint size = filesize;
		const char* c = &bytes[0];

		GLuint s = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(s,1, &c, &size);
		glCompileShader(s);
		glAttachShader(p,s);
		m_printInfoLogs(p, s);
	}
	std::cerr << "done." << std::endl;

	glLinkProgram(p);

	return p;
}


GLuint m_CreateFragmentShader(const char* shader_path)
{
	// create a fragment shader (the vertex shader is using the fixed-function pipeline)
	GLuint p = glCreateProgram();
	std::cerr << "compiling " << shader_path << "...";
	{
		std::ifstream ifs(shader_path, std::ios::binary);
		ifs.seekg(0, std::ios::end);
		std::ifstream::pos_type filesize = ifs.tellg();
		ifs.seekg(0, std::ios::beg);
		std::vector<char> bytes(filesize);
		ifs.read(&bytes[0], filesize);
		GLint size = filesize;
		const char* c = &bytes[0];

		GLuint s = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(s,1, &c, &size);
		glCompileShader(s);
		glAttachShader(p,s);
		m_printInfoLogs(p, s);
	}
	std::cerr << "done." << std::endl;

	glLinkProgram(p);

	return p;
}


void m_idle(void)
{
	glutPostRedisplay();
}


GLuint m_CreateTexture(const int TextureIndex, GLenum format, const int BufferSize)
{
	GLuint Texture;

	glActiveTexture(GL_TEXTURE0 + TextureIndex);
	glGenTextures(1, &Texture);
	glBindTexture(GL_TEXTURE_2D, Texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexImage2D(GL_TEXTURE_2D, 0, format, BufferSize, BufferSize, 0, GL_LUMINANCE, GL_FLOAT, 0);

	return Texture;
}


int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - ImageResolution) / 2, (glutGet(GLUT_SCREEN_HEIGHT) - ImageResolution) / 2);
	glutInitWindowSize(ImageResolution, ImageResolution);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH);
	glutCreateWindow(argv[0]);
	glutDisplayFunc(m_display);
	glutIdleFunc(m_idle);

	// create shaders
	PSDraw = m_CreateFragmentShader("draw.fs");
	PSHash = m_CreateFragmentShader("hash.fs");
	
	//tigra: для всех фотонов проверить на видимость(или еще для чего-то QueryIntersection
	PSProgressiveUpdate = m_CreateFragmentShader("progressive.fs");
	
	//tigra: для пикселов всех разделить flux на радиус с пи и т.д.
	PSRadianceEstimate = m_CreateFragmentShader("re.fs");
	PVSScatter = m_CreateFullShader("scatter.vs", "scatter.fs");
	PVSCorrection = m_CreateFullShader("correction.vs", "correction.fs");
	PSEyeRayTrace = m_CreateFragmentShader("eyeraytrace.fs");
	PSPhotonTrace = m_CreateFragmentShader("photontrace.fs");
	PSMax = m_CreateFragmentShader("max.fs");
	PSMin = m_CreateFragmentShader("min.fs");
	PSSum = m_CreateFragmentShader("sum.fs");

	// create textures
	PhotonHashTexture = m_CreateTexture(0, GL_RGBA32F_ARB, HashResolution);
	PhotonCorrectionTexture = m_CreateTexture(1, GL_RGBA32F_ARB, HashResolution);

	QueryNormalTexture = m_CreateTexture(2, GL_RGBA32F_ARB, ImageResolution);
	QueryPositionTexture = m_CreateTexture(3, GL_RGBA32F_ARB, ImageResolution);
	RandomEyeRayTexture = m_CreateTexture(4, GL_RGBA32F_ARB, ImageResolution);
	RandomPhotonTexture = m_CreateTexture(4, GL_RGBA32F_ARB, PhotonBufferSize);
	QueryEmissionPhotonCountTexture = m_CreateTexture(5, GL_RGBA32F_ARB, ImageResolution);
	QueryFluxRadiusTexture = m_CreateTexture(6, GL_RGBA32F_ARB, ImageResolution);
	QueryReflectanceTexture = m_CreateTexture(7, GL_RGBA32F_ARB, ImageResolution);

	PhotonIndexTexture = m_CreateTexture(8, GL_RGBA32F_ARB, PhotonBufferSize);
	PhotonFluxTexture = m_CreateTexture(9, GL_RGBA32F_ARB, PhotonBufferSize);
	PhotonPositionTexture = m_CreateTexture(10, GL_RGBA32F_ARB, PhotonBufferSize);
	PhotonDirectionTexture = m_CreateTexture(11, GL_RGBA32F_ARB, PhotonBufferSize);

	PhotonIntersectionTexture = m_CreateTexture(15, GL_RGBA32F_ARB, PhotonBufferSize);
	PhotonEmittedFlagTexture = m_CreateTexture(15, GL_RGBA32F_ARB, PhotonBufferSize);
	QueryIntersectionTexture = m_CreateTexture(15, GL_RGBA32F_ARB, ImageResolution);

	// buffer for computing min/max/average
	MinMaxAveTextureQuery = m_CreateTexture(12, GL_RGBA32F_ARB, ImageResolution);
	glGenFramebuffersEXT(1, &MinMaxAveSurfaceQuery);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, MinMaxAveSurfaceQuery);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, MinMaxAveTextureQuery, 0);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	MinMaxAveTexturePhoton = m_CreateTexture(12, GL_RGBA32F_ARB, PhotonBufferSize);
	glGenFramebuffersEXT(1, &MinMaxAveSurfacePhoton);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, MinMaxAveSurfacePhoton);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, MinMaxAveTexturePhoton, 0);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	// create FBOs
	// precomputed hash values
	glGenFramebuffersEXT(1, &PhotonIndexSurface);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, PhotonIndexSurface);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, PhotonIndexTexture, 0);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	// hash buffer
	glGenFramebuffersEXT(1, &PhotonHashSurface);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, PhotonHashSurface);
	glGenRenderbuffersEXT(1, &PhotonHashDepthBuffer);
	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, PhotonHashDepthBuffer);
	glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, HashResolution, HashResolution);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, PhotonHashTexture, 0);
	glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, PhotonHashDepthBuffer);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	// hash-count buffer
	glGenFramebuffersEXT(1, &PhotonCorrectionSurface);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, PhotonCorrectionSurface);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, PhotonCorrectionTexture, 0);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	// eye ray intersection data
	GLenum EyeRayTraceBuffers[] = {GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT, GL_COLOR_ATTACHMENT2_EXT, GL_COLOR_ATTACHMENT3_EXT, GL_COLOR_ATTACHMENT4_EXT, GL_COLOR_ATTACHMENT5_EXT};
	glGenFramebuffersEXT(1, &EyeRayTraceSurface);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, EyeRayTraceSurface);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, QueryPositionTexture, 0);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT, GL_TEXTURE_2D, QueryReflectanceTexture, 0);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT2_EXT, GL_TEXTURE_2D, QueryNormalTexture, 0);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT3_EXT, GL_TEXTURE_2D, RandomEyeRayTexture, 0);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT4_EXT, GL_TEXTURE_2D, QueryIntersectionTexture, 0);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT5_EXT, GL_TEXTURE_2D, QueryEmissionPhotonCountTexture, 0);
	glDrawBuffers(6, EyeRayTraceBuffers);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	// photon data
	GLenum PhotonRayTraceBuffers[] = {GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT, GL_COLOR_ATTACHMENT2_EXT, GL_COLOR_ATTACHMENT3_EXT, GL_COLOR_ATTACHMENT4_EXT, GL_COLOR_ATTACHMENT5_EXT};
	glGenFramebuffersEXT(1, &PhotonRayTraceSurface);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, PhotonRayTraceSurface);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, PhotonPositionTexture, 0);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT, GL_TEXTURE_2D, PhotonFluxTexture, 0);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT2_EXT, GL_TEXTURE_2D, PhotonDirectionTexture, 0);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT3_EXT, GL_TEXTURE_2D, RandomPhotonTexture, 0);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT4_EXT, GL_TEXTURE_2D, PhotonIntersectionTexture, 0);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT5_EXT, GL_TEXTURE_2D, PhotonEmittedFlagTexture, 0);
	glDrawBuffers(6, PhotonRayTraceBuffers);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	// measurement points
	GLenum QueryBuffers[] = {GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT};
	glGenFramebuffersEXT(1, &QueryPointSurface);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, QueryPointSurface);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, QueryFluxRadiusTexture, 0);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT, GL_TEXTURE_2D, QueryEmissionPhotonCountTexture, 0);
	glDrawBuffers(2, QueryBuffers);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	// create a VBO
	glGenBuffersARB(1, &FragmentsVBO);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, FragmentsVBO);
	std::vector<float> VBOData;
	for (int j = 0; j < PhotonBufferSize; j++)
	{
		for (int i = 0; i < PhotonBufferSize; i++)
		{
			VBOData.push_back(2.0f * ((i + 0.5f) / float(PhotonBufferSize)) - 1.0f);
			VBOData.push_back(2.0f * ((j + 0.5f) / float(PhotonBufferSize)) - 1.0f);
		}
	}
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(float) * 2 * PhotonBufferSize * PhotonBufferSize, (void*)&VBOData[0], GL_STATIC_DRAW_ARB);

	// initialize misc data
	FrameCount = 0;
	glClearColor(0.0, 0.0, 0.0, 0.0);
	
	CanonicalCameraPosition = TVector3(0.0f, 0.0f, 13.0f);
	
	// load mesh data
	if(argc>1)
	{		
		gpurt.mesh.LoadOBJ(argv[1], TVector3(0.0f, 0.0f, 0.0f), 1);
		CanonicalCameraPosition = TVector3(0.0f, 0.0f, 5.0f);
	}
	else
	gpurt.mesh.LoadOBJ("cornell_metal.obj", TVector3(0.0f, 0.0f, 0.0f), 0.01f);

	FieldOfView = 45.0f;
	LookAtPosition = TVector3(0.0f, 0.0f, 0.0f);
	gpurt.camera.Set(CanonicalCameraPosition, LookAtPosition, ImageResolution, ImageResolution, FieldOfView);


	// precalcuation (BVH construction) for mesh
	std::cerr << "building BVH...";
	gpurt.PrecalculateMeshData();
	std::cerr << "done" << std::endl;

	if (gpurt.mesh.lightsCDF.size() == 0)
	{
		std::cerr << "no light source is defined, use constant illumination" << std::endl;
	}

	// enter the main loop
	std::cerr << "start rendering..." << std::endl;
	RandomizeTextures();
	glClear(GL_COLOR_BUFFER_BIT);
	startedTime = glutGet(GLUT_ELAPSED_TIME);
	glutMainLoop();

	// release the resources of gpurt
	gpurt.Release();

	// delete things
	glDeleteBuffersARB(1, &FragmentsVBO);

	glDeleteProgram(PSDraw);
	glDeleteProgram(PSHash);
	glDeleteProgram(PSProgressiveUpdate);
	glDeleteProgram(PSRadianceEstimate);
	glDeleteProgram(PVSScatter);
	glDeleteProgram(PVSCorrection);
	glDeleteProgram(PSEyeRayTrace);
	glDeleteProgram(PSPhotonTrace);

	glDeleteProgram(PSMax);
	glDeleteProgram(PSMin);
	glDeleteProgram(PSSum);

	glDeleteTextures(1, &QueryPositionTexture);
	glDeleteTextures(1, &QueryNormalTexture);
	glDeleteTextures(1, &QueryEmissionPhotonCountTexture);
	glDeleteTextures(1, &QueryFluxRadiusTexture);
	glDeleteTextures(1, &QueryReflectanceTexture);
	glDeleteTextures(1, &QueryIntersectionTexture);

	glDeleteTextures(1, &PhotonIndexTexture);
	glDeleteTextures(1, &PhotonFluxTexture);
	glDeleteTextures(1, &PhotonPositionTexture);
	glDeleteTextures(1, &PhotonDirectionTexture);
	glDeleteTextures(1, &RandomPhotonTexture);
	glDeleteTextures(1, &RandomEyeRayTexture);
	glDeleteTextures(1, &PhotonHashTexture);
	glDeleteTextures(1, &PhotonCorrectionTexture);
	glDeleteTextures(1, &PhotonIntersectionTexture);
	glDeleteTextures(1, &PhotonEmittedFlagTexture);

	glDeleteTextures(1, &MinMaxAveTextureQuery);
	glDeleteFramebuffersEXT(1, &MinMaxAveSurfaceQuery);
	glDeleteTextures(1, &MinMaxAveTexturePhoton);
	glDeleteFramebuffersEXT(1, &MinMaxAveSurfacePhoton);

	glDeleteFramebuffersEXT(1, &PhotonHashSurface);
	glDeleteFramebuffersEXT(1, &PhotonCorrectionSurface);
	glDeleteFramebuffersEXT(1, &PhotonIndexSurface);
	glDeleteFramebuffersEXT(1, &QueryPointSurface);
	glDeleteFramebuffersEXT(1, &EyeRayTraceSurface);
	glDeleteFramebuffersEXT(1, &PhotonRayTraceSurface);
	glDeleteRenderbuffersEXT(1, &PhotonHashDepthBuffer);

	return 0;
}
