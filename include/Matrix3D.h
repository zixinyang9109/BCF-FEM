// Matrix3D.h : Header declaration for the 3D matrix transform class
//

#pragma once
#include <fstream>

//#include "CommonDefs.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////

class Matrix3D
{
	////////////////////////////////////////////
	// Member Variables
	public:
		float Data[4][4];


	////////////////////////////////////////////
	// Member Functions
	public:
		Matrix3D();
		~Matrix3D();

		friend ostream& operator<<(ostream& os, const Matrix3D& Src);
		float* operator[](int Row);

		Matrix3D& operator=(const Matrix3D& Src);
		Matrix3D operator*(const Matrix3D& Src) const;
		Matrix3D operator+(Matrix3D& Src) const;
		Matrix3D operator-(Matrix3D& Src) const;

		Matrix3D& operator=(float Val);
		Matrix3D operator*(float Val) const;
		Matrix3D operator/(float Val) const;
		Matrix3D operator+(float Val) const;
		Matrix3D operator-(float Val) const;

		static Matrix3D GetIdentity(void);
		static Matrix3D GetTrans(float Tx, float Ty, float Tz);
		static Matrix3D GetdTrans(float Tx, float Ty, float Tz);
		static Matrix3D GetScale(float Sx, float Sy, float Sz);

		static Matrix3D GetRotX(float Angle);
		static Matrix3D GetRotY(float Angle);
		static Matrix3D GetRotZ(float Angle);

		static Matrix3D GetdRotX(float Angle, float sign);
		static Matrix3D GetdRotY(float Angle, float sign);
		static Matrix3D GetdRotZ(float Angle, float sign);

		static Matrix3D GetTranspose(Matrix3D& Src);
		static Matrix3D GetInverse(Matrix3D& Src);
		static Matrix3D GetInverseRT(Matrix3D& Src);

		static Matrix3D PerspectZ(float Dz);
		static Matrix3D PerspectZ(float Dx, float Dy, float Dz);
		static Matrix3D PerspectX(float Dx);
		static Matrix3D Perspect(float v[4], float n[4]);
		void VectorMult(float Src[4], float Dst[4]) const;
		static void Homogenize(float Src[4], float Dst[4]);
		void Proj2D(float Src[3], float Dst[2]) const;
		void RotationMatrixToEuler(float &angleX, float &angleY, float &angleZ);
};

