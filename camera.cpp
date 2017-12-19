#include "camera.h"


void TCamera::Set(const TVector3& _origin, const TVector3& _lookat, const int _width, const int _height, const float _fov)
{
	this->origin = _origin;
	this->lookat = _lookat;

	this->width = _width;
	this->height = _height;
	this->distance = float(_height) / (2.0f * tanf((_fov / 2.0f) * (3.141592f / 180.0f)));

	const TVector3 tv = TVector3(0.0f, 1.0f, 0.0f);

	this->w = (this->lookat - this->origin); this->w = this->w.normalize();
	this->u = tv % this->w; this->u = this->u.normalize();
	this->v = this->w % this->u; this->v = this->v.normalize();
}
