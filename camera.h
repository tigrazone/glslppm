#ifndef _CAMERA_H
#define	_CAMERA_H


#include "vector.h"
#include <math.h>


class TCamera
{
	public:
		TVector3 origin, lookat;

		int width, height;
		float distance;

		TVector3 u, v, w;

		void Set(const TVector3& _origin, const TVector3& _lookat, const int _width, const int _height, const float _fov);
	private:
};


#endif
