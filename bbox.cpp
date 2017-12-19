#include "bbox.h"


int TBBox::LongestAxis() const {
	const TVector3 length = this->max - this->min;

	if ((length.x > length.y) && (length.x > length.z)) {
		return 0;
	}
	else if (length.y > length.z) {
		return 1;
	}
	else {
		return 2;
	}
}


float TBBox::Area() const {
	const TVector3 length = this->max - this->min;
	return 2.0f * (length.x * length.y + length.y * length.z + length.z * length.x);
}


void TBBox::Expand(const TVector3& p) {
	this->max = this->max.Maximize(p);
	this->min = this->min.Minimize(p);
}


void TBBox::Initialize() {
	this->min = TVector3(1e+38f, 1e+38f, 1e+38f);
	this->max = TVector3(-1e+38f, -1e+38f, -1e+38f);
}
