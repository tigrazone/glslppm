#ifndef _BBOX_H
#define _BBOX_H


#include "vector.h"


class TBBox {
public:
	TVector3 min, max;

	float Area() const;
	void Expand(const TVector3 &p);
	void Initialize();
	int LongestAxis() const;
};


#endif
