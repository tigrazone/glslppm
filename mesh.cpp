#include "mesh.h"


void get_indices(char *word, int *vindex, int *tindex, int *nindex)
{
	char *null = " ";
	char *ptr;
	char *tp;
	char *np;

	tp = null;
	np = null;

	for (ptr = word; *ptr != '\0'; ptr++) 
	{
		if (*ptr == '/') {
			if (tp == null)
			{
				tp = ptr + 1;
			}
			else
			{
				np = ptr + 1;
			}

			*ptr = '\0';
		}
	}

	*vindex = atoi(word);
	*tindex = atoi(tp);
	*nindex = atoi(np);
}


std::vector<TMaterial> mtls;
void LoadPPM(const char *fname, const int i)
{
	FILE *fp;

	fp = fopen(fname, "rb");

	fscanf(fp, "P6\n");
	fscanf(fp, "%d %d\n", &mtls[i].textureWidth, &mtls[i].textureHeight);
	fscanf(fp, "255\n");
	mtls[i].texture.resize(mtls[i].textureWidth * mtls[i].textureHeight * 3);
	size_t s = fread(&mtls[i].texture[0], 1, mtls[i].textureWidth * mtls[i].textureHeight * 3, fp);
	if (s != mtls[i].texture.size())
	{
		return;
	}
	fclose(fp);
}


void LoadMTL(const std::string fileName)
{
	FILE * fp = fopen(fileName.c_str(),"r");

	TMaterial mtl;
	mtl.texture.clear();
	char line[81];
	while (fgets(line, 80, fp) != NULL)
	{
		float r, g, b, s;
		std::string lineStr;
		lineStr = line;
		int i = mtls.size();

		if (lineStr.compare(0, 6, "newmtl", 0, 6) == 0)
		{
			lineStr.erase(0, 7);
			mtl.name = lineStr;
			mtl.isTextured = false;
		}
		else if (lineStr.compare(0, 2, "Ka", 0, 2) == 0)
		{
			lineStr.erase(0, 3);
			sscanf(lineStr.c_str(), "%f %f %f\n", &r, &g, &b);
			mtl.Ka = TVector3(r, g, b);
		}
		else if (lineStr.compare(0, 2, "Kd", 0, 2) == 0)
		{
			lineStr.erase(0, 3);
			sscanf(lineStr.c_str(), "%f %f %f\n", &r, &g, &b);
			mtl.Kd = TVector3(r, g, b);
		}
		else if (lineStr.compare(0, 2, "Ks", 0, 2) == 0)
		{
			lineStr.erase(0, 3);
			sscanf(lineStr.c_str(), "%f %f %f\n", &r, &g, &b);
			mtl.Ks = TVector3(r, g, b);
		}
		else if (lineStr.compare(0, 2, "Ns", 0, 2) == 0)
		{
			lineStr.erase(0, 3);
			sscanf(lineStr.c_str(), "%f\n", &s);
			mtl.Ns = s;
			mtls.push_back(mtl);
			mtls[i].texture.clear();
		}
		else if (lineStr.compare(0, 6, "map_Kd", 0, 6) == 0)
		{
			lineStr.erase(0, 7);
			lineStr.erase(lineStr.size() - 1, 1);
			mtls[i - 1].isTextured = true;
			LoadPPM(lineStr.c_str(), i - 1);
		}
	}

	fclose(fp);
}


void ParseOBJ(char* fileName, int &nVertices, float **vertices, float **normals, float **texcoords, int &nIndices, int **indices, int **materials)
{
	FILE * fp = fopen(fileName,"r");
	int nv = 0, nn = 0, nf = 0, nt = 0;
	char line[81];

	while (fgets( line, 80, fp ) != NULL)
	{
		std::string lineStr;
		lineStr = line;

		if (lineStr.compare(0, 6, "mtllib", 0, 6) == 0)
		{
			lineStr.erase(0, 7);
			lineStr.erase(lineStr.size() - 1, 1);
			LoadMTL(lineStr);
		}

		if (line[0] == 'v')
		{
			if (line[1] == 'n')
			{
				nn++;
			}
			else if (line[1] == 't')
			{
				nt++;
			}
			else
			{
				nv++;
			}
		}
		else if (line[0] == 'f')
		{
			nf++;
		}
	}
	fseek(fp, 0, 0);

	float *n = new float[3 * (nn > nf ? nn : nf)];
	float *v = new float[3 * nv];
	float *t = new float[2 * nt];

	int *vInd = new int[3 * nf];
	int *nInd = new int[3 * nf];
	int *tInd = new int[3 * nf];
	int *mInd = new int[nf];

	int nvertices = 0;
	int nnormals = 0;
	int ntexcoords = 0;
	int nindices  = 0;
	int ntriangles = 0;
	bool noNormals = false;
	bool noTexCoords = false;
	bool noMaterials = true;
	int cmaterial = 0;

	while (fgets( line, 80, fp ) != NULL)
	{
		std::string lineStr;
		lineStr = line;

		if (line[0] == 'v')
		{
			if (line[1] == 'n')
			{
				float x, y, z;
				sscanf(&line[2], "%f %f %f\n", &x, &y, &z);
				float l = sqrt(x * x + y * y + z * z);
				x = x / l;
				y = y / l;
				z = z / l;
				n[nnormals] = x;
				nnormals++;
				n[nnormals] = y;
				nnormals++;
				n[nnormals] = z;
				nnormals++;
			}
			else if (line[1] == 't')
			{
				float u, v;
				sscanf(&line[2], "%f %f\n", &u, &v);
				t[ntexcoords] = u;
				ntexcoords++;
				t[ntexcoords] = v;
				ntexcoords++;
			}
			else
			{
				float x, y, z;
				sscanf( &line[1], "%f %f %f\n", &x, &y, &z);
				v[nvertices] = x;
				nvertices++;
				v[nvertices] = y;
				nvertices++;
				v[nvertices] = z;
				nvertices++;
			}
		}
		if (lineStr.compare(0, 6, "usemtl", 0, 6) == 0)
		{
			lineStr.erase(0, 7);

			if (mtls.size() != 0)
			{
				for (unsigned int i = 0; i < mtls.size(); i++)
				{
					if (lineStr.compare(mtls[i].name) == 0)
					{
						cmaterial = i;
						noMaterials = false;
						break;
					}
				}
			}

		}
		else if (line[0] == 'f') 
		{
			char s1[32], s2[32], s3[32];
			int vI, tI, nI;
			sscanf(&line[1], "%s %s %s\n", s1, s2, s3);

			mInd[ntriangles] = cmaterial;

			// indices for first vertex
			get_indices(s1, &vI, &tI, &nI);
			vInd[nindices] = vI - 1;
			if (nI)
			{
				nInd[nindices] = nI - 1;
			}
			else
			{
				noNormals = true;
			}

			if (tI)
			{
				tInd[nindices] = tI - 1;
			}
			else
			{
				noTexCoords = true;
			}
			nindices++;

			// indices for second vertex
			get_indices(s2, &vI, &tI, &nI);
			vInd[nindices] = vI - 1;
			if (nI)
			{
				nInd[nindices] = nI - 1;
			}
			else
			{
				noNormals = true;
			}

			if (tI)
			{
				tInd[nindices] = tI - 1;
			}
			else
			{
				noTexCoords = true;
			}
			nindices++;

			// indices for third vertex
			get_indices(s3, &vI, &tI, &nI);
			vInd[nindices] = vI - 1;
			if (nI)
			{
				nInd[nindices] = nI - 1;
			}
			else
			{
				noNormals = true;
			}

			if (tI)
			{
				tInd[nindices] = tI - 1;
			}
			else
			{
				noTexCoords = true;
			}
			nindices++;

			ntriangles++;
		}
	}

	// we don't support separate indices for normals, vertices, and texture coordinates. 
	*vertices = new float[ntriangles*9];
	if (!noNormals) 
	{
		*normals = new float[ntriangles*9];
	}
	else 
	{
		*normals = 0;
	}

	if (!noTexCoords)
	{
		*texcoords = new float[ntriangles*6];
	}
	else
	{
		*texcoords = 0;
	}

	if (!noMaterials)
	{
		*materials = new int[ntriangles];
	}
	else
	{
		*materials = 0;
	}

	*indices = new int[ntriangles*3];
	nVertices = ntriangles*3;
	nIndices = ntriangles*3;

	for (int i = 0; i < ntriangles; i++)
	{
		if (!noMaterials)
		{
			(*materials)[i] = mInd[i];
		}

		(*indices)[3 * i] = 3 * i;
		(*indices)[3 * i + 1] = 3 * i + 1;
		(*indices)[3 * i + 2] = 3 * i + 2;

		(*vertices)[9 * i] = v[3 * vInd[3 * i]];
		(*vertices)[9 * i + 1] = v[3 * vInd[3 * i] + 1];
		(*vertices)[9 * i + 2] = v[3 * vInd[3 * i] + 2];

		(*vertices)[9 * i + 3] = v[3 * vInd[3 * i + 1]];
		(*vertices)[9 * i + 4] = v[3 * vInd[3 * i + 1] + 1];
		(*vertices)[9 * i + 5] = v[3 * vInd[3 * i + 1] + 2];

		(*vertices)[9 * i + 6] = v[3 * vInd[3 * i + 2]];
		(*vertices)[9 * i + 7] = v[3 * vInd[3 * i + 2] + 1];
		(*vertices)[9 * i + 8] = v[3 * vInd[3 * i + 2] + 2];

		if(!noNormals)
		{
			(*normals)[9 * i] = n[3 * nInd[3 * i]];
			(*normals)[9 * i + 1] = n[3 * nInd[3 * i]+1];
			(*normals)[9 * i + 2] = n[3 * nInd[3 * i]+2];

			(*normals)[9 * i + 3] = n[3*nInd[3 * i + 1]];
			(*normals)[9 * i + 4] = n[3*nInd[3 * i + 1] + 1];
			(*normals)[9 * i + 5] = n[3*nInd[3 * i + 1] + 2];

			(*normals)[9 * i + 6] = n[3*nInd[3 * i + 2]];
			(*normals)[9 * i + 7] = n[3*nInd[3 * i + 2] + 1];
			(*normals)[9 * i + 8] = n[3*nInd[3 * i + 2] + 2];
		}

		if(!noTexCoords)
		{
			(*texcoords)[6 * i ] = t[2*tInd[3 * i]];
			(*texcoords)[6 * i + 1] = t[2*tInd[3 * i] + 1];

			(*texcoords)[6 * i + 2] = t[2*tInd[3 * i + 1]];
			(*texcoords)[6 * i + 3] = t[2*tInd[3 * i + 1] + 1];

			(*texcoords)[6 * i + 4] = t[2*tInd[3 * i + 2]];
			(*texcoords)[6 * i + 5] = t[2*tInd[3 * i + 2] + 1];
		}

	}

	fclose(fp);

	delete[] n;
	delete[] v;
	delete[] t;
	delete[] nInd;
	delete[] vInd;
	delete[] tInd;
	delete[] mInd;
}


void TMesh::Release()
{
	this->triangles.clear();
	this->materials.clear();
}


void TMesh::CalculateBBox()
{
	this->bbox.Initialize();
	for (unsigned int i = 0; i < this->triangles.size(); i++)
	{
		this->triangles[i].bbox.Initialize();
		this->triangles[i].bbox.Expand(this->triangles[i].positions[0]);
		this->triangles[i].bbox.Expand(this->triangles[i].positions[1]);
		this->triangles[i].bbox.Expand(this->triangles[i].positions[2]);

		this->triangles[i].centroid = (this->triangles[i].positions[0] + this->triangles[i].positions[1] + this->triangles[i].positions[2]) * (1.0f / 3.0f);

		this->bbox.Expand(this->triangles[i].positions[0]);
		this->bbox.Expand(this->triangles[i].positions[1]);
		this->bbox.Expand(this->triangles[i].positions[2]);
	}
}


void TMesh::PrepareLightSources()
{
	this->lightsArea = 0.0f;
	this->lightsCDF.clear();
	this->lightsIndices.clear();

	for (unsigned int i = 0; i < this->triangles.size(); i++)
	{
		if (this->materials[this->triangles[i].idMaterial].brdf == -1)
		{
			const TVector3 Edge0 = this->triangles[i].positions[1] - this->triangles[i].positions[0];
			const TVector3 Edge1 = this->triangles[i].positions[2] - this->triangles[i].positions[0];
			const float a = 0.5f * (Edge0 % Edge1).length();
			this->lightsCDF.push_back(a);
			this->lightsIndices.push_back(i);
			this->lightsArea += a;
		}
	}

	if (this->lightsCDF.size() == 0) return;

	for (unsigned int i = 1; i < this->lightsCDF.size(); i++)
	{
		this->lightsCDF[i] = this->lightsCDF[i] + this->lightsCDF[i - 1];
	}

	for (unsigned int i = 0; i < this->lightsCDF.size(); i++)
	{
		this->lightsCDF[i] = this->lightsCDF[i] / this->lightsArea;
	}
}


void TMesh::LoadOBJ(char* FileName, TVector3 Position, float Scale)
{
	int nVertices;
	float *vertices;
	float *normals;
	float *texcoords;
	int nIndices;
	int *indices;
	int *matid;
	bool haveLightSource = false;

	ParseOBJ(FileName, nVertices, &vertices, &normals, &texcoords, nIndices, &indices, &matid);

	this->triangles.resize(nIndices / 3);

	if (matid != NULL)
	{
		for (unsigned int i = 0; i < mtls.size(); i++)
		{
			this->materials.push_back(mtls[i]);
			if (mtls[i].isTextured)
			{
				this->materials[i].texture.clear();
				for (unsigned int j = 0; j < mtls[i].texture.size(); j++)
				{
					this->materials[i].texture.push_back(mtls[i].texture[j]);
				}
			}
		}
		for (unsigned int i = 0; i < mtls.size(); i++)
		{
			// Lambertian
			this->materials[i].brdf = 0;
			this->materials[i].eta = 1.7f;
			this->materials[i].specularity = 1.0f;

			if (mtls[i].Ns == 100.0f)
			{
				if (mtls[i].Ks.dot(mtls[i].Ks) == 3.0f)
				{
					// mirror
					this->materials[i].brdf = 1;
				}
				else if (mtls[i].Ks.dot(mtls[i].Ks) > 0.0f)
				{
					// plastic
					this->materials[i].brdf = 3;
				}
			}
			else
			{
				this->materials[i].specularity = std::max(mtls[i].Ns / 100.0f, 0.5f);
				if (mtls[i].Ks.dot(mtls[i].Ks) == 3.0f)
				{
					// glossy mirror
					this->materials[i].brdf = 4;
				}
				else if (mtls[i].Ks.dot(mtls[i].Ks) > 0.0f)
				{
					// glossy plastic
					this->materials[i].brdf = 6;
				}
			}

			if (mtls[i].name.compare(0, 5, "glass", 0, 5) == 0)
			{
			this->materials[i].eta = mtls[i].Ks.dot(mtls[i].Ks) / 3.0f + 1.0f;
				if (mtls[i].Ns == 100.0f)
				{
					// glass
					this->materials[i].brdf = 2;
				}
				else 
				{
					// glossy glass
					this->materials[i].specularity = std::max(mtls[i].Ns / 100.0f, 0.5f);
					this->materials[i].brdf = 5;
				}
			}

			if (mtls[i].Ka.dot(mtls[i].Ka) > 0.0f)
			{
				// light source
				this->materials[i].brdf = -1;
				haveLightSource = true;
			}

			if ((this->materials[i].brdf == 0) || (this->materials[i].brdf == 3))
			{
				mtls[i].Kd = mtls[i].Kd * 0.9f;
			}

			if ((this->materials[i].brdf == 2) || (this->materials[i].brdf == 5))
			{
				mtls[i].Kd.x = sqrtf(mtls[i].Kd.x);
				mtls[i].Kd.y = sqrtf(mtls[i].Kd.y);
				mtls[i].Kd.z = sqrtf(mtls[i].Kd.z);
			}

			this->materials[i].color = mtls[i].Kd;

			if (mtls[i].name.compare(0, 3, "sss", 0, 3) == 0)
			{
				this->materials[i].brdf = 7;
			}
		}
	}
	else
	{
		// use default
		TMaterial mtl;
		mtl.isTextured = false;
		mtl.brdf = 0;
		mtl.Kd = TVector3(0.75f, 0.75f, 0.75f);
		this->materials.push_back(mtl);
	}

	for (unsigned int i = 0; i < this->triangles.size(); i++)
	{
		const int v0 = indices[i * 3 + 0];
		const int v1 = indices[i * 3 + 1];
		const int v2 = indices[i * 3 + 2];

		this->triangles[i].positions[0] = TVector3(vertices[v0 * 3 + 0], vertices[v0 * 3 + 1], vertices[v0 * 3 + 2]);
		this->triangles[i].positions[1] = TVector3(vertices[v1 * 3 + 0], vertices[v1 * 3 + 1], vertices[v1 * 3 + 2]);
		this->triangles[i].positions[2] = TVector3(vertices[v2 * 3 + 0], vertices[v2 * 3 + 1], vertices[v2 * 3 + 2]);
		this->triangles[i].positions[0] = this->triangles[i].positions[0] * Scale + Position;
		this->triangles[i].positions[1] = this->triangles[i].positions[1] * Scale + Position;
		this->triangles[i].positions[2] = this->triangles[i].positions[2] * Scale + Position;

		if (normals != NULL)
		{
			this->triangles[i].normals[0] = TVector3(normals[v0 * 3 + 0], normals[v0 * 3 + 1], normals[v0 * 3 + 2]);
			this->triangles[i].normals[1] = TVector3(normals[v1 * 3 + 0], normals[v1 * 3 + 1], normals[v1 * 3 + 2]);
			this->triangles[i].normals[2] = TVector3(normals[v2 * 3 + 0], normals[v2 * 3 + 1], normals[v2 * 3 + 2]);
		}
		else
		{
			// no normal data, calculate the normal for a polygon
			const TVector3 e0 = this->triangles[i].positions[1] - this->triangles[i].positions[0];
			const TVector3 e1 = this->triangles[i].positions[2] - this->triangles[i].positions[0];
			const TVector3 n = (e0 % e1).normalize();

			this->triangles[i].normals[0] = n;
			this->triangles[i].normals[1] = n;
			this->triangles[i].normals[2] = n;
		}

		// material id
		this->triangles[i].idMaterial = 0;
		if (matid != NULL)
		{
			// read texture coordinates
			if ((texcoords != NULL) && mtls[matid[i]].isTextured)
			{
				this->triangles[i].texcoords[0] = TVector2(texcoords[v0 * 2 + 0], texcoords[v0 * 2 + 1]);
				this->triangles[i].texcoords[1] = TVector2(texcoords[v1 * 2 + 0], texcoords[v1 * 2 + 1]);
				this->triangles[i].texcoords[2] = TVector2(texcoords[v2 * 2 + 0], texcoords[v2 * 2 + 1]);
			}
			else
			{
				this->triangles[i].texcoords[0] = TVector2(1.0e+30f, 1.0e+30f);
				this->triangles[i].texcoords[1] = TVector2(1.0e+30f, 1.0e+30f);
				this->triangles[i].texcoords[2] = TVector2(1.0e+30f, 1.0e+30f);
			}

			this->triangles[i].idMaterial = matid[i];
		}
		else
		{
			this->triangles[i].texcoords[0] = TVector2(1.0e+30f, 1.0e+30f);
			this->triangles[i].texcoords[1] = TVector2(1.0e+30f, 1.0e+30f);
			this->triangles[i].texcoords[2] = TVector2(1.0e+30f, 1.0e+30f);
		}
	}

	this->CalculateBBox();
	if (haveLightSource) this->PrepareLightSources();

	delete[] vertices;
	delete[] normals;
	delete[] texcoords;
	delete[] indices;
	delete[] matid;
}
