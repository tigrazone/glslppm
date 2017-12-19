#include "gpurt.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>


void CGPURT::PrecalculateMeshData()
{
	glDeleteTextures(1, &this->CubeTextureBBoxRootIndices);
	glDeleteTextures(1, &this->textureBVHs);
	glDeleteTextures(1, &this->textureTriangles);

	this->BuildBVH();
	this->CreateTextures();
}


GLuint CGPURT::CreateTexture(const int TextureIndex, GLenum format, const int Width, const int Height) 
{
	GLuint Texture;

	glActiveTexture(GL_TEXTURE0 + TextureIndex);
	glGenTextures(1, &Texture);
	glBindTexture(GL_TEXTURE_2D, Texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexImage2D(GL_TEXTURE_2D, 0, format, Width, Height, 0, GL_LUMINANCE, GL_FLOAT, 0);

	return Texture;
}
 
float sgn(const float v)
{
	if (v > 0.0f)
	{
		return 1.0f;
	}
	else
	{
		return -1.0f;
	}
}

void CGPURT::CreateTextures()
{
	int k = 0;

	const int numTriangles = this->mesh.triangles.size();
	std::vector<TVector2> triangleTexcoord(numTriangles);

	int PolyDataSizeX = int(sqrtf(numTriangles * 6.0f) + 1.0f);
	int PolyDataSizeY = PolyDataSizeX;

	this->PolygonDataStride.x = 1.0f / float(PolyDataSizeX);
	this->PolygonDataStride.y = 1.0f / float(PolyDataSizeY);


	// --------------- creating textures for positions ---------------
	TVector4* PolyTexSys = new TVector4[PolyDataSizeX * PolyDataSizeY];
	this->textureTriangles = CreateTexture(0, GL_RGBA32F_ARB, PolyDataSizeX, PolyDataSizeY);
	{
		k = 0;
		while (k < (numTriangles * 6))
		{
			int p = k / 6;
			int i = k % PolyDataSizeX;
			int j = k / PolyDataSizeX;
			triangleTexcoord[p].x = ((float)i + 0.5f) / (float)PolyDataSizeX;
			triangleTexcoord[p].y = ((float)j + 0.5f) / (float)PolyDataSizeY;

			for (int v = 0; v <= 2; v++)
			{
				PolyTexSys[k].x = this->mesh.triangles[p].positions[v].x;
				PolyTexSys[k].y = this->mesh.triangles[p].positions[v].y;
				PolyTexSys[k].z = this->mesh.triangles[p].positions[v].z;
				PolyTexSys[k].w = this->mesh.triangles[p].texcoords[v].x;
				k++;
			}

			for (int v = 0; v <= 2; v++)
			{
				PolyTexSys[k].x = this->mesh.triangles[p].normals[v].x;
				PolyTexSys[k].y = this->mesh.triangles[p].normals[v].y;
				PolyTexSys[k].z = sgn(this->mesh.triangles[p].normals[v].z) * (this->mesh.triangles[p].idMaterial + 1);
				PolyTexSys[k].w = this->mesh.triangles[p].texcoords[v].y;
				k++;
			}

		}
	}
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, PolyDataSizeX, PolyDataSizeY, 0, GL_RGBA, GL_FLOAT, PolyTexSys);
	delete[] PolyTexSys;



	// --------------- material data ---------------
	if (this->mesh.materials.size() == 0)
	{
		// no material data found, use the default
		this->mesh.materials.resize(1);
		this->mesh.materials[0].color = TVector3(0.5f);
		this->mesh.materials[0].brdf = 0;
		this->mesh.materials[0].eta = 1.0f;
	}

	int MaterialDataNumVec4s = 2;
	// (vec3(r,g,b), 0.0)
	// (BRDF, eta, vec2(0.0))
	this->MaterialDataStride = 1.0f / (float)this->mesh.materials.size();
	TVector4* MaterialTexSys = new TVector4[this->mesh.materials.size() * MaterialDataNumVec4s];
	this->TextureMaterials = CreateTexture(6, GL_RGBA32F_ARB, this->mesh.materials.size() * MaterialDataNumVec4s, 1);
		k = 0;
		for (unsigned int i = 0; i <= this->mesh.materials.size() - 1; i++)
		{
			MaterialTexSys[k].x = this->mesh.materials[i].color.x;
			MaterialTexSys[k].y = this->mesh.materials[i].color.y;
			MaterialTexSys[k].z = this->mesh.materials[i].color.z;
			MaterialTexSys[k].w = this->mesh.materials[i].eta;
			k++;

			MaterialTexSys[k].x = (float)this->mesh.materials[i].brdf;
			MaterialTexSys[k].y = this->mesh.materials[i].specularity;
			MaterialTexSys[k].z = 0.0f;
			MaterialTexSys[k].w = 0.0f;
			k++;
		}
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, this->mesh.materials.size() * MaterialDataNumVec4s, 1, 0, GL_RGBA, GL_FLOAT, MaterialTexSys);
	delete[] MaterialTexSys;






	// --------------- light source data ---------------
	if (this->mesh.lightsCDF.size() > 0)
	{
		TVector4* LightSourcesTexSys = new TVector4[this->mesh.lightsCDF.size()];
		this->TextureLightSources = CreateTexture(6, GL_RGBA32F_ARB, this->mesh.lightsCDF.size(), 1);
			k = 0;
			for (unsigned int i = 0; i <= this->mesh.lightsCDF.size() - 1; i++)
			{
				LightSourcesTexSys[k].x = triangleTexcoord[this->mesh.lightsIndices[i]].x;
				LightSourcesTexSys[k].y = triangleTexcoord[this->mesh.lightsIndices[i]].y;
				LightSourcesTexSys[k].z = this->mesh.lightsCDF[i];
				LightSourcesTexSys[k].w = 0.0f;
				k++;
			}
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, this->mesh.lightsCDF.size(), 1, 0, GL_RGBA, GL_FLOAT, LightSourcesTexSys);
		delete[] LightSourcesTexSys;
	}


	glEnable(GL_TEXTURE_3D);
	glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, &this->VolumeTextureTextures);
	glBindTexture(GL_TEXTURE_3D, this->VolumeTextureTextures);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);


	// resampling textures into fixed resolution to simplify the implementation (not absolutely necessary)
	const int TexRes = 512;
	std::vector<unsigned char> VolumeTextureTexturesSys(TexRes * TexRes * this->mesh.materials.size() * 3);

	for (unsigned int z = 0; z < this->mesh.materials.size(); z++)
	{
		for (unsigned int y = 0; y < TexRes; y++)
		{
			for (unsigned int x = 0; x < TexRes; x++)
			{
				const int dst = 3 * (x + y * TexRes + z * TexRes * TexRes);
				const int src = 3 * ((x * this->mesh.materials[z].textureWidth / TexRes) + (y * this->mesh.materials[z].textureHeight / TexRes) * this->mesh.materials[z].textureWidth);
				if (this->mesh.materials[z].isTextured)
				{ 
					VolumeTextureTexturesSys[dst + 0] = this->mesh.materials[z].texture[src + 0];
					VolumeTextureTexturesSys[dst + 1] = this->mesh.materials[z].texture[src + 1];
					VolumeTextureTexturesSys[dst + 2] = this->mesh.materials[z].texture[src + 2];
				}
				else 
				{
					VolumeTextureTexturesSys[dst + 0] = 255;
					VolumeTextureTexturesSys[dst + 1] = 0;
					VolumeTextureTexturesSys[dst + 2] = 255;
				}
			}
		}
	}

	glTexImage3D(GL_TEXTURE_3D, 0, GL_RGB8, TexRes, TexRes, this->mesh.materials.size(), 0, GL_RGB, GL_UNSIGNED_BYTE, &VolumeTextureTexturesSys[0]);
	VolumeTextureTexturesSys.clear();


	// --------------- recording texture coordinates for each bounding box data ---------------
	int BBoxNum = this->bvh.nodesNum;
	TBBoxGPU* BBox[6];
	TVector2* BBoxTexCoord[6];

	for (int kk = 0; kk <= 5; kk++)
	{
		BBox[kk] = new TBBoxGPU[BBoxNum];
		BBoxTexCoord[kk] = new TVector2[BBoxNum];
	}


	BBoxDataSizeX = int(sqrtf(float(BBoxNum)) + 1.0f);
	BBoxDataSizeY = BBoxDataSizeX;

	// root selector
	glEnable(GL_TEXTURE_CUBE_MAP);
	glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, &this->CubeTextureBBoxRootIndices);
	glBindTexture(GL_TEXTURE_CUBE_MAP, this->CubeTextureBBoxRootIndices);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	for (int i = 0; i < 6; i++) 
	{ 
		TVector4 VData;
		if (i == 0)
		{
			VData.x = (((float)this->BBoxDataSizeX * 0.0f + 0.5f) / ((float)this->BBoxDataSizeX * 2.0f));
			VData.y = (((float)this->BBoxDataSizeY * 0.0f + 0.5f) / ((float)this->BBoxDataSizeY * 4.0f));
			VData.z = (-1.0f);
			VData.w = (-1.0f);
		}
		else if (i == 1)
		{
			VData.x = (((float)this->BBoxDataSizeX * 1.0f + 0.5f) / ((float)this->BBoxDataSizeX * 2.0f));
			VData.y = (((float)this->BBoxDataSizeY * 0.0f + 0.5f) / ((float)this->BBoxDataSizeY * 4.0f));
			VData.z = (-1.0f);
			VData.w = (-1.0f);
		}
		else if (i == 2)
		{
			VData.x = (((float)this->BBoxDataSizeX * 0.0f + 0.5f) / ((float)this->BBoxDataSizeX * 2.0f));
			VData.y = (((float)this->BBoxDataSizeY * 1.0f + 0.5f) / ((float)this->BBoxDataSizeY * 4.0f));
			VData.z = (-1.0f);
			VData.w = (-1.0f);
		}
		else if (i == 3)
		{
			VData.x = (((float)this->BBoxDataSizeX * 1.0f + 0.5f) / ((float)this->BBoxDataSizeX * 2.0f));
			VData.y = (((float)this->BBoxDataSizeY * 1.0f + 0.5f) / ((float)this->BBoxDataSizeY * 4.0f));
			VData.z = (-1.0f);
			VData.w = (-1.0f);
		}
		else if (i == 4)
		{
			VData.x = (((float)this->BBoxDataSizeX * 0.0f + 0.5f) / ((float)this->BBoxDataSizeX * 2.0f));
			VData.y = (((float)this->BBoxDataSizeY * 2.0f + 0.5f) / ((float)this->BBoxDataSizeY * 4.0f));
			VData.z = (-1.0f);
			VData.w = (-1.0f);
		}
		else
		{
			VData.x = (((float)this->BBoxDataSizeX * 1.0f + 0.5f) / ((float)this->BBoxDataSizeX * 2.0f));
			VData.y = (((float)this->BBoxDataSizeY * 2.0f + 0.5f) / ((float)this->BBoxDataSizeY * 4.0f));
			VData.z = (-1.0f);
			VData.w = (-1.0f);
		}

		glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGBA32F_ARB, 1, 1, 0, GL_RGBA, GL_FLOAT, &(VData.x));
	}



	// --------------- record min, max, hit, miss and index of triangle ---------------
	// recording 2D texture coordinates
	int cOffSetX, cOffSetY;
	for (int kk = 0; kk <= 5; kk++)
	{
		if (kk == 0) 
		{
			cOffSetX = this->BBoxDataSizeX * 0;
			cOffSetY = this->BBoxDataSizeY * 0;
		}
		else if (kk == 1)
		{
			cOffSetX = this->BBoxDataSizeX * 1;
			cOffSetY = this->BBoxDataSizeY * 0;
		}
		else if (kk == 2)
		{
			cOffSetX = this->BBoxDataSizeX * 0;
			cOffSetY = this->BBoxDataSizeY * 1;
		}
		else if (kk == 3)
		{
			cOffSetX = this->BBoxDataSizeX * 1;
			cOffSetY = this->BBoxDataSizeY * 1;
		}
		else if (kk == 4) 
		{
			cOffSetX = this->BBoxDataSizeX * 0;
			cOffSetY = this->BBoxDataSizeY * 2;
		}
		else
		{
			cOffSetX = this->BBoxDataSizeX * 1;
			cOffSetY = this->BBoxDataSizeY * 2;
		}

		// generate bbox texture coordinates
		k = 0;
		for (int j = 0; j <= this->BBoxDataSizeY - 1; j++)
		{
			for (int i = 0; i <= this->BBoxDataSizeX - 1; i++)
			{
				if (k < BBoxNum) 
				{
					BBoxTexCoord[kk][k].x = (float)(i + cOffSetX + 0.5f) / ((float)this->BBoxDataSizeX * 2.0f);
					BBoxTexCoord[kk][k].y = (float)(j + cOffSetY + 0.5f) / ((float)this->BBoxDataSizeY * 4.0f);
				}
				k++;
			}
		}

		for (int i = 0; i <= BBoxNum - 1; i++)
		{
			// --------------- the canonical BVH ---------------
			BBox[kk][i].min = this->bvh.nodes[kk][i].bbox.min;
			BBox[kk][i].max = this->bvh.nodes[kk][i].bbox.max;

			if ((i + 1) >= BBoxNum)
			{
				BBox[kk][i].hit = TVector2(-1.0f, -1.0f);
			}
			else
			{
				BBox[kk][i].hit = BBoxTexCoord[kk][i + 1];
			}

			BBox[kk][i].tri = -1;
			if (this->bvh.nodes[kk][i].isLeaf)
			{
				if (this->bvh.nodes[kk][i].idTriangle != -1)
				{
					BBox[kk][i].tri = this->bvh.nodes[kk][i].idTriangle;
				}
				else
				{
					BBox[kk][i].tri = -1;
				}
			}

			if ((this->bvh.nodes[kk][i].idMiss == -1) || (this->bvh.nodes[kk][i].idMiss >= BBoxNum))
			{
				BBox[kk][i].miss = TVector2(-1.0f, -1.0f);
			}
			else
			{
				BBox[kk][i].miss = BBoxTexCoord[kk][this->bvh.nodes[kk][i].idMiss];
			}

			// --------------- other BVHs ---------------
			if (kk != 0) 
			{
				BBox[kk][i].min = this->bvh.nodes[0][i].bbox.min;
				BBox[kk][i].max = this->bvh.nodes[0][i].bbox.max;

				if ((i + 1) >= BBoxNum)
				{
					BBox[kk][i].hit = TVector2(-1.0f, -1.0f);
				}
				else
				{
					BBox[kk][i].hit = BBoxTexCoord[kk][this->bvh.nodes[kk][i + 1].idBase];
				}

				BBox[kk][i].tri = -1;
				if (this->bvh.nodes[0][i].isLeaf)
				{
					if (this->bvh.nodes[0][i].idTriangle != -1)
					{
						BBox[kk][i].tri = this->bvh.nodes[0][i].idTriangle;
					}
					else
					{
						BBox[kk][i].tri = -1;
					}
				}

				if ((this->bvh.nodes[kk][i].idMiss == -1) || (this->bvh.nodes[kk][i].idMiss >= BBoxNum))
				{
					BBox[kk][i].miss = TVector2(-1.0f, -1.0f);
				}
				else
				{
					BBox[kk][i].miss = BBoxTexCoord[kk][this->bvh.nodes[kk][this->bvh.nodes[kk][i].idMiss].idBase];
				}
			}

		}
	}
	// --------------------------------------------------




	// --------------- write textures for BVH (min and max, hit and miss) ---------------
	// note that the size of the texture is (BBoxDataSizeX * 2, BBoxDataSizeY * 4)
	TVector4* TextureBVHSys = new TVector4[(this->BBoxDataSizeX * 2) * (this->BBoxDataSizeY * 4)];
	this->textureBVHs = CreateTexture(7, GL_RGBA32F_ARB, this->BBoxDataSizeX * 2, this->BBoxDataSizeY * 4);
	{
		int* temp_bbox_index = new int[BBoxNum];

		// --------------- for all BVHs ---------------
		for (int kk = 0; kk <= 5; kk++)
		{
			// calculate index to the canonical BVH data
			for (int i = 0; i <= BBoxNum - 1; i++)
			{
				temp_bbox_index[this->bvh.nodes[kk][i].idBase] = i;
			}

			// --------------- add offset for each BVH data ---------------
			k = 0;
			int bbhitmiss = 0;


			if (kk == 0)
			{
				bbhitmiss = 0;
			}
			else if (kk == 1)
			{
				bbhitmiss = this->BBoxDataSizeX;
			}
			else if (kk == 2)
			{
				bbhitmiss = (this->BBoxDataSizeX * 2) * this->BBoxDataSizeY;
			}
			else if (kk == 3)
			{
				bbhitmiss = (this->BBoxDataSizeX * 2) * this->BBoxDataSizeY + this->BBoxDataSizeX;
			}
			else if (kk == 4)
			{
				bbhitmiss = (this->BBoxDataSizeX * 2) * this->BBoxDataSizeY * 2;
			}
			else if (kk == 5)
			{
				bbhitmiss = (this->BBoxDataSizeX * 2) * this->BBoxDataSizeY * 2 + this->BBoxDataSizeX;
			}

			// --------------- write hit and miss links ---------------
			for (int j = 0; j <= this->BBoxDataSizeY - 1; j++)
			{
				for (int i = 0; i <= this->BBoxDataSizeX - 1; i++)
				{
					if (k < BBoxNum)
					{
						TextureBVHSys[bbhitmiss].x = (BBox[kk][temp_bbox_index[k]].hit.x);
						TextureBVHSys[bbhitmiss].y = (BBox[kk][temp_bbox_index[k]].hit.y);
						TextureBVHSys[bbhitmiss].z = (BBox[kk][temp_bbox_index[k]].miss.y);
						TextureBVHSys[bbhitmiss].w = (BBox[kk][temp_bbox_index[k]].miss.x);
					}

					bbhitmiss++;
					k++;
				}
				// skip to the correct next scanline
				bbhitmiss = bbhitmiss + this->BBoxDataSizeX;
			}
		}


		// --------------- write bounding box min and x coordinates of triangle index ---------------
		k = 0;
		int bbhitmiss = (this->BBoxDataSizeX * 2) * this->BBoxDataSizeY * 3;

		for (int j = 0; j <= this->BBoxDataSizeY - 1; j++)
		{
			for (int i = 0; i <= this->BBoxDataSizeX - 1; i++)
			{
				TextureBVHSys[bbhitmiss] = TVector4(-1, -1, -1, 0);
				if (k < BBoxNum)
				{
					TextureBVHSys[bbhitmiss].x = (BBox[0][k].min.x);
					TextureBVHSys[bbhitmiss].y = (BBox[0][k].min.y);
					TextureBVHSys[bbhitmiss].z = (BBox[0][k].min.z);
					if (BBox[0][k].tri != -1)
					{
						TextureBVHSys[bbhitmiss].w = (triangleTexcoord[BBox[0][k].tri].x);
					}
					else 
					{
						TextureBVHSys[bbhitmiss].w = (-1.0f);
					}
				}
				bbhitmiss++;
				k++;
			}
			bbhitmiss = bbhitmiss + this->BBoxDataSizeX * 1;
		}


		// --------------- write bounding box max and y coordinates of triangle index ---------------
		k = 0;
		bbhitmiss = (this->BBoxDataSizeX * 2) * this->BBoxDataSizeY * 3 + this->BBoxDataSizeX;

		for (int j = 0; j <= this->BBoxDataSizeY - 1; j++)
		{
			for (int i = 0; i <= this->BBoxDataSizeX - 1; i++)
			{
				TextureBVHSys[bbhitmiss] = TVector4(1, 1, 1, 0);
				if (k < BBoxNum)
				{
					TextureBVHSys[bbhitmiss].x = (BBox[0][k].max.x);
					TextureBVHSys[bbhitmiss].y = (BBox[0][k].max.y);
					TextureBVHSys[bbhitmiss].z = (BBox[0][k].max.z);
					if (BBox[0][k].tri != -1)
					{
						TextureBVHSys[bbhitmiss].w = (triangleTexcoord[BBox[0][k].tri].y);
					}
					else 
					{
						TextureBVHSys[bbhitmiss].w = (-1.0f);
					}
				}
				bbhitmiss++;
				k++;
			}
			bbhitmiss = bbhitmiss + this->BBoxDataSizeX * 1;
		}
		delete[] temp_bbox_index;
	}
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, this->BBoxDataSizeX * 2, this->BBoxDataSizeY * 4, 0, GL_RGBA, GL_FLOAT, TextureBVHSys);
	delete[] TextureBVHSys;

	delete[] BBox[0];
	delete[] BBox[1];
	delete[] BBox[2];
	delete[] BBox[3];
	delete[] BBox[4];
	delete[] BBox[5];

	delete[] BBoxTexCoord[0];
	delete[] BBoxTexCoord[1];
	delete[] BBoxTexCoord[2];
	delete[] BBoxTexCoord[3];
	delete[] BBoxTexCoord[4];
	delete[] BBoxTexCoord[5];
}


void CGPURT::BuildBVH()
{
	this->bvh.Build(this->mesh);
}


void CGPURT::Release()
{
	glDeleteTextures(1, &this->textureBVHs);
	glDeleteTextures(1, &this->VolumeTextureTextures);

	this->mesh.Release();

	glDeleteTextures(1, &this->textureTriangles);
	glDeleteTextures(1, &this->CubeTextureBBoxRootIndices);
	glDeleteTextures(1, &this->TextureMaterials);
	glDeleteTextures(1, &this->TextureLightSources);
}