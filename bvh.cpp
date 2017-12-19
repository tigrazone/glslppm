#include "bvh.h"


static int tnodeNum;


void sortAxis(TMesh& mesh, int* obj_index, char axis, int li, int ri)
{
	int i = li;
	int j = ri;

	const float pivot = mesh.triangles[obj_index[(li + ri) / 2]].centroid[axis];
	for (;;)
	{
		while (mesh.triangles[obj_index[i]].centroid[axis] < pivot) i++;
		while (mesh.triangles[obj_index[j]].centroid[axis] > pivot) j--;
		if (i >= j) break;

		const int temp = obj_index[i];
		obj_index[i] = obj_index[j];
		obj_index[j] = temp;

		i++;
		j--;
	}

	if (li < (i - 1)) sortAxis(mesh, obj_index, axis, li, i - 1);
	if ((j + 1) <  ri) sortAxis(mesh, obj_index, axis, j + 1, ri);
}


int splitBVH(TMesh& mesh, int* obj_index, int obj_num, TBBox& bbox, int Level, int face, TBVHNode** mnode)
{
	int* obj_index_L;
	int* obj_index_R;

	// --------------- leaf node ---------------
	// subdivision is done until we have one node - this is to simplify the implementation on a GPU, but can be changed of course.
	if (obj_num <= 1) 
	{
		tnodeNum++;
		const int temp_id = tnodeNum - 1;

		mnode[face][temp_id].bbox = bbox;
		mnode[face][temp_id].isLeaf = true;

		if (obj_num != 0)
		{
			mnode[face][temp_id].idTriangle = obj_index[0];
		}
		else
		{
			mnode[face][temp_id].idTriangle = -1;
		}

		return temp_id;
	}


	// --------------- internal node ---------------
	int bestAxis = 0;
	int bestIndex = 0;
	float bestCost = 1e+30f;
	TBBox bestBBoxL, bestBBoxR;

	// obvious case
	if (obj_num == 2)
	{
		// divide the node into two nodes
		obj_index_L = new int[1];
		obj_index_R = new int[1];

		obj_index_L[0] = obj_index[0];
		obj_index_R[0] = obj_index[1];

		bestBBoxL = mesh.triangles[obj_index[0]].bbox;
		bestBBoxR = mesh.triangles[obj_index[1]].bbox;
	}
	else
	{
		// --------------- use exact SAH (surface area heuristic) by sorting ---------------
		int* sorted_obj_index = new int[obj_num];
		float* leftArea = new float[obj_num];
		TBBox* lbbox = new TBBox[obj_num];

		// for all axes
		for (int k = 0; k <= 2; k++)
		{
			sortAxis(mesh, obj_index, k, 0, obj_num - 1);

			// calculate area of bounding boxes for left sweeping
			TBBox bboxL, bboxR;
			bboxL.Initialize();
			for (int i = 0; i <= obj_num - 1; i++)
			{
				const int ii = obj_index[i];
				bboxL.Expand(mesh.triangles[ii].bbox.min);
				bboxL.Expand(mesh.triangles[ii].bbox.max);
				leftArea[i] = bboxL.Area();
				lbbox[i] = bboxL;
			}

			// calculate SAH by right sweeping
			int triNum = obj_num - 1;
			bboxR.Initialize();
			for (int j = (obj_num - 2); j >= 0; j--)
			{
				const int ii = obj_index[j + 1];

				bboxR.Expand(mesh.triangles[ii].bbox.min);
				bboxR.Expand(mesh.triangles[ii].bbox.max);

				const float tempCost = (float)triNum * leftArea[j] + (float)(obj_num - triNum) * bboxR.Area();
				if (tempCost < bestCost) 
				{
					bestCost = tempCost;
					bestAxis = k;
					bestIndex = j;
					bestBBoxL = lbbox[bestIndex];
					bestBBoxR = bboxR;
				}

				triNum--;
			}

			// use the best axis
			if (bestAxis == k)
			{
				for (int i = 0; i <= obj_num - 1; i++)
				{
					sorted_obj_index[i] = obj_index[i];
				}
			}
		}

		// divide the node into two nodes
		obj_index_L = new int[bestIndex + 1];
		obj_index_R = new int[obj_num - (bestIndex + 1)];

		for (int i = 0; i <= bestIndex; i++)
		{
			obj_index_L[i] = sorted_obj_index[i];
		}
		for (int i = bestIndex + 1; i <= obj_num - 1; i++)
		{
			obj_index_R[i - (bestIndex + 1)] = sorted_obj_index[i];
		}

		delete[] sorted_obj_index;
	}

	// it is not a leaf node
	tnodeNum++;
	const int temp_id = tnodeNum - 1;
	mnode[face][temp_id].bbox = bbox;
	mnode[face][temp_id].isLeaf = false;

	// follow canonical condition to make BVH
	if (bestBBoxL.min.x < bestBBoxR.min.x)
	{
		mnode[face][temp_id].idLeft = splitBVH(mesh, obj_index_L, bestIndex + 1, bestBBoxL, Level + 1, face, mnode);
		mnode[face][temp_id].idRight = splitBVH(mesh, obj_index_R, obj_num - (bestIndex + 1), bestBBoxR, Level + 1, face, mnode);
	}
	else
	{
		mnode[face][temp_id].idLeft = splitBVH(mesh, obj_index_R, obj_num - (bestIndex + 1), bestBBoxR, Level + 1, face, mnode);
		mnode[face][temp_id].idRight = splitBVH(mesh, obj_index_L, bestIndex + 1, bestBBoxL, Level + 1, face, mnode);
	}

	delete[] obj_index_L;
	delete[] obj_index_R;

	return temp_id;
}



void ReorderNodes(TMesh& mesh, int face, int index, TBVHNode** mnode)
{
	if (index < 0) return;
	if ((unsigned int)tnodeNum == (mesh.triangles.size() * 2)) return;

	tnodeNum++;
	int temp_id = tnodeNum - 1;
	mnode[face][temp_id] = mnode[6][index];
	mnode[face][temp_id].idBase = index;

	if (mnode[6][index].isLeaf) return;

	ReorderNodes(mesh, face, mnode[6][index].idLeft, mnode);
	ReorderNodes(mesh, face, mnode[6][index].idRight, mnode);
}


int ReorderTree(TMesh& mesh, int face, int index, TBVHNode** mnode)
{
	if (mnode[6][index].isLeaf)
	{
		tnodeNum++;
		return tnodeNum - 1;
	}

	tnodeNum++;
	int temp_id = tnodeNum - 1;
	mnode[face][temp_id].idLeft = ReorderTree(mesh, face, mnode[6][index].idLeft, mnode);
	mnode[face][temp_id].idRight = ReorderTree(mesh, face, mnode[6][index].idRight, mnode);
	return temp_id;
}


void SetLeftMissLinks(int id, int idParent, int face, TBVHNode** mnode)
{
	if (mnode[face][id].isLeaf) 
	{
		mnode[face][id].idMiss = id + 1;
		return;
	}

	mnode[face][id].idMiss = mnode[face][idParent].idRight;

	SetLeftMissLinks(mnode[face][id].idLeft, id, face, mnode);
	SetLeftMissLinks(mnode[face][id].idRight, id, face, mnode);
}


void SetRightMissLinks(int id, int idParent, int face, TBVHNode** mnode)
{
	if (mnode[face][id].isLeaf)
	{
		mnode[face][id].idMiss = id + 1;
		return; 
	}

	if (mnode[face][idParent].idRight == id)
	{
		mnode[face][id].idMiss = mnode[face][idParent].idMiss;
	}

	SetRightMissLinks(mnode[face][id].idLeft, id, face, mnode);
	SetRightMissLinks(mnode[face][id].idRight, id, face, mnode);
}


void TBVH::Build(TMesh& mesh)
{
	this->nodes = new TBVHNode*[7];

	// build six BVHs for all the canonical directions
	for (int face = 0; face <= 5; face++)
	{
		// initialize nodes
		const int obj_num = mesh.triangles.size();
		int* obj_index = new int[obj_num];
		for (int i = 0; i <= obj_num - 1; i++)
		{
			obj_index[i] = i;
		}
		tnodeNum = 0;
		this->nodes[face] = new TBVHNode[obj_num * 2];
		for (int i = 0; i <= obj_num * 2 - 1; i++)
		{
			this->nodes[face][i].idMiss = -1;
			this->nodes[face][i].idBase = i;
		}

		if (face == 0)
		{
			// canonical BVH (optimal BVH for rays go to positive-x directions)
			splitBVH(mesh, obj_index, obj_num, mesh.bbox, 0, face, this->nodes);
			this->nodesNum = tnodeNum;

			// initialize temporary BVH nodes
			this->nodes[6] = new TBVHNode[obj_num * 2];
			for (int i = 0; i <= obj_num * 2 - 1; i++)
			{
				this->nodes[6][i].idMiss = -1;
			}
		}
		else
		{
			// other BVHs
			for (int i = 0; i <= this->nodesNum - 1; i++)
			{
				this->nodes[6][i] = this->nodes[0][i];
			}

			// swap indices if certain conditions are met for each BVH
			for (int i = 0; i <= this->nodesNum - 1; i++)
			{
				if (this->nodes[6][i].isLeaf) continue;

				if ((face == 1) && (this->nodes[6][this->nodes[6][i].idLeft].bbox.max.x > this->nodes[6][this->nodes[6][i].idRight].bbox.max.x)) continue;
				if ((face == 2) && (this->nodes[6][this->nodes[6][i].idLeft].bbox.min.y < this->nodes[6][this->nodes[6][i].idRight].bbox.min.y)) continue;
				if ((face == 3) && (this->nodes[6][this->nodes[6][i].idLeft].bbox.max.y > this->nodes[6][this->nodes[6][i].idRight].bbox.max.y)) continue;
				if ((face == 4) && (this->nodes[6][this->nodes[6][i].idLeft].bbox.min.z < this->nodes[6][this->nodes[6][i].idRight].bbox.min.z)) continue;
				if ((face == 5) && (this->nodes[6][this->nodes[6][i].idLeft].bbox.max.z > this->nodes[6][this->nodes[6][i].idRight].bbox.max.z)) continue;

				const int temp = this->nodes[6][i].idLeft;
				this->nodes[6][i].idLeft = this->nodes[6][i].idRight;
				this->nodes[6][i].idRight = temp;
			}

			// rebuilding BVH
			tnodeNum = 0;
			ReorderNodes(mesh, face, 0, this->nodes);
			tnodeNum = 0;
			ReorderTree(mesh, face, 0, this->nodes);
		}

		// threading BVH (making miss links) 
		this->nodes[face][0].idMiss = -1;
		SetLeftMissLinks(0, 0, face, this->nodes);
		this->nodes[face][0].idMiss = -1;
		SetRightMissLinks(0, 0, face, this->nodes);
		this->nodes[face][0].idMiss = -1;
	}
}
