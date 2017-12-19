#ifndef _BVH_H
#define _BVH_H


#include "vector.h"
#include "mesh.h"
#include "bbox.h"


class TBVHNode
{
	public:
		TBBox bbox;

		bool isLeaf;

		int idLeft, idRight, idTriangle;
		int idMiss, idBase;

	private:
};


class TBVH
{
	public:
		// nodes
		TBVHNode** nodes;
		int nodesNum;

		void Build(TMesh& mesh);
	private:
};


#endif
