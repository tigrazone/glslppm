#ifndef _VECTOR_H
#define _VECTOR_H


#include <math.h>


struct TVector2
{
	float x, y;
	TVector2(float x_ = 0.0f, float y_ = 0.0f) {x = x_; y = y_;}
};


struct TVector4
{
	float x, y, z, w;
	TVector4(float x_ = 0.0f, float y_ = 0.0f, float z_ = 0.0f, float w_ = 0.0f) {x = x_; y = y_; z = z_; w = w_;}
	inline float& operator[](const int i)       {if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; else return w;}
	inline float  operator[](const int i) const {if (i == 0) return x; else if (i == 1) return y; else if (i == 2) return z; else return w;}
};


struct TVector3
{
	float x, y, z;
	TVector3(float x_ = 0.0f, float y_ = 0.0f, float z_ = 0.0f) {x = x_; y = y_; z = z_;}
	inline TVector3 operator+(const TVector3 &b) const {return TVector3(x+b.x, y+b.y, z+b.z);}
	inline TVector3 operator-(const TVector3 &b) const {return TVector3(x-b.x, y-b.y, z-b.z);}
	inline TVector3 operator+(const float b) const { return TVector3(x + b, y + b, z + b); }
	inline TVector3 operator-(const float b) const { return TVector3(x - b, y - b, z - b); }
	inline TVector3 operator*(const float b) const { return TVector3(x * b, y * b, z * b); }
	inline TVector3 operator/(const float b) const { return TVector3(x / b, y / b, z / b); }
	inline TVector3 mul(const TVector3 &b) const {return TVector3(x * b.x, y * b.y , z * b.z);}
	inline TVector3 const normalize() {return (*this) * (1.0f / sqrt(x*x+y*y+z*z));}
	inline float dot(const TVector3 &b) const {return x * b.x + y * b.y + z * b.z;}
	inline float length() const {return sqrtf(this->dot(*this));}
	inline TVector3 operator%(const TVector3&b) const {return TVector3(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
	inline TVector3 Maximize(const TVector3 &b) const {return TVector3(x>b.x?x:b.x, y>b.y?y:b.y, z>b.z?z:b.z);}
	inline TVector3 Minimize(const TVector3 &b) const {return TVector3(x<b.x?x:b.x, y<b.y?y:b.y, z<b.z?z:b.z);}
	inline float& operator[](const int i)       {if (i == 0) return x; else if (i == 1) return y; else return z;}
	inline float  operator[](const int i) const {if (i == 0) return x; else if (i == 1) return y; else return z;}
};


#endif
