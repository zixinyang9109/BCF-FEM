// Matrix3D.cpp : Function definitions for the 3D matrix transform class
//

#include <math.h>
#include "Matrix3D.h"

///////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////

Matrix3D::Matrix3D()
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			Data[i][j] = 0.0f;
		}
	}
}

///////////////////////////////////////////////////////////////////////////
// Destructor.
///////////////////////////////////////////////////////////////////////////

Matrix3D::~Matrix3D()
{
}

///////////////////////////////////////////////////////////////////////////
// Matrix output
///////////////////////////////////////////////////////////////////////////

ostream& operator<<(ostream& os, const Matrix3D& Src)
{
	return os << Src.Data[0][0] << " " << Src.Data[0][1] << " " << Src.Data[0][2] << " " << Src.Data[0][3] << endl
	          << Src.Data[1][0] << " " << Src.Data[1][1] << " " << Src.Data[1][2] << " " << Src.Data[1][3] << endl
	          << Src.Data[2][0] << " " << Src.Data[2][1] << " " << Src.Data[2][2] << " " << Src.Data[2][3] << endl
	          << Src.Data[3][0] << " " << Src.Data[3][1] << " " << Src.Data[3][2] << " " << Src.Data[3][3] << endl << endl;
}

///////////////////////////////////////////////////////////////////////////
// Matrix assignment to another matrix
///////////////////////////////////////////////////////////////////////////

float* Matrix3D::operator[](int Row)
{
	return Data[Row];
}

///////////////////////////////////////////////////////////////////////////
// Matrix assignment to another matrix
///////////////////////////////////////////////////////////////////////////

Matrix3D& Matrix3D::operator=(const Matrix3D& Src)
{
	for (int j=0; j<4; j++)
		for (int i=0; i<4; i++)
			Data[j][i] = Src.Data[j][i];

	return *this;
}

///////////////////////////////////////////////////////////////////////////
// Matrix multiplication
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::operator*(const Matrix3D& Src) const
{
Matrix3D tmp;
	
	for (int j=0; j<4; j++)
		for (int i=0; i<4; i++)
			tmp.Data[j][i] = Data[j][0]*Src.Data[0][i] + Data[j][1]*Src.Data[1][i] + Data[j][2]*Src.Data[2][i] + Data[j][3]*Src.Data[3][i];

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Matrix addition
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::operator+(Matrix3D& Src) const
{
Matrix3D tmp;
	
	for (int j=0; j<4; j++)
		for (int i=0; i<4; i++)
			tmp.Data[j][i] = Data[j][i] + Src.Data[j][i];

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Matrix subtraction
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::operator-(Matrix3D& Src) const
{
Matrix3D tmp;
	
	for (int j=0; j<4; j++)
		for (int i=0; i<4; i++)
			tmp.Data[j][i] = Data[j][i] - Src.Data[j][i];

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Matrix assignment to an immediate value
///////////////////////////////////////////////////////////////////////////

Matrix3D& Matrix3D::operator=(float Val)
{
	for (int j=0; j<4; j++)
		for (int i=0; i<4; i++)
				Data[j][i] = Val;

	return *this;
}

///////////////////////////////////////////////////////////////////////////
// Matrix multiplication by an immediate value
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::operator*(float Val) const
{
Matrix3D tmp;
	
	for (int j=0; j<4; j++)
		for (int i=0; i<4; i++)
			tmp.Data[j][i] = Data[j][i] * Val;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Matrix division by an immediate value
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::operator/(float Val) const
{
Matrix3D tmp;
	
	for (int j=0; j<4; j++)
		for (int i=0; i<4; i++)
			tmp.Data[j][i] = Data[j][i] / Val;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Matrix addition with an immediate value
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::operator+(float Val) const
{
Matrix3D tmp;
	
	for (int j=0; j<4; j++)
		for (int i=0; i<4; i++)
			tmp.Data[j][i] = Data[j][i] + Val;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Matrix subtraction with an immediate value
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::operator-(float Val) const
{
Matrix3D tmp;
	
	for (int j=0; j<4; j++)
		for (int i=0; i<4; i++)
			tmp.Data[j][i] = Data[j][i] - Val;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Set to an identity matrix.
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::GetIdentity(void)
{
Matrix3D tmp;

	tmp.Data[0][0] = 1; tmp.Data[0][1] = 0; tmp.Data[0][2] = 0; tmp.Data[0][3] = 0;
	tmp.Data[1][0] = 0; tmp.Data[1][1] = 1; tmp.Data[1][2] = 0; tmp.Data[1][3] = 0;
	tmp.Data[2][0] = 0; tmp.Data[2][1] = 0; tmp.Data[2][2] = 1; tmp.Data[2][3] = 0;
	tmp.Data[3][0] = 0; tmp.Data[3][1] = 0; tmp.Data[3][2] = 0; tmp.Data[3][3] = 1;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Set to a translation matrix.
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::GetTrans(float Tx, float Ty, float Tz)
{
Matrix3D tmp;

	tmp.Data[0][0] = 1; tmp.Data[0][1] = 0; tmp.Data[0][2] = 0; tmp.Data[0][3] = Tx;
	tmp.Data[1][0] = 0; tmp.Data[1][1] = 1; tmp.Data[1][2] = 0; tmp.Data[1][3] = Ty;
	tmp.Data[2][0] = 0; tmp.Data[2][1] = 0; tmp.Data[2][2] = 1; tmp.Data[2][3] = Tz;
	tmp.Data[3][0] = 0; tmp.Data[3][1] = 0; tmp.Data[3][2] = 0; tmp.Data[3][3] = 1;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Set to a translation matrix.
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::GetdTrans(float Tx, float Ty, float Tz)
{
	Matrix3D tmp;

	tmp.Data[0][0] = 0; tmp.Data[0][1] = 0; tmp.Data[0][2] = 0; tmp.Data[0][3] = Tx;
	tmp.Data[1][0] = 0; tmp.Data[1][1] = 0; tmp.Data[1][2] = 0; tmp.Data[1][3] = Ty;
	tmp.Data[2][0] = 0; tmp.Data[2][1] = 0; tmp.Data[2][2] = 0; tmp.Data[2][3] = Tz;
	tmp.Data[3][0] = 0; tmp.Data[3][1] = 0; tmp.Data[3][2] = 0; tmp.Data[3][3] = 0;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Set to a scaling matrix.
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::GetScale(float Sx, float Sy, float Sz)
{
Matrix3D tmp;

	tmp.Data[0][0] = Sx; tmp.Data[0][1] = 0;  tmp.Data[0][2] = 0;  tmp.Data[0][3] = 0;
	tmp.Data[1][0] = 0;  tmp.Data[1][1] = Sy; tmp.Data[1][2] = 0;  tmp.Data[1][3] = 0;
	tmp.Data[2][0] = 0;  tmp.Data[2][1] = 0;  tmp.Data[2][2] = Sz; tmp.Data[2][3] = 0;
	tmp.Data[3][0] = 0;  tmp.Data[3][1] = 0;  tmp.Data[3][2] = 0;  tmp.Data[3][3] = 1;

	return tmp;
}

/////////////////////////////////////////////////////////////////////////////////////
//The following rotation matrices are vector rotation matrices. They corrspond to the 
//rotation of a vector(object)relative to fixed axes.  
////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// Set to a matrix to rotate around the x axis.
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::GetRotX(float Angle)
{
Matrix3D tmp;

	float ca = cos(Angle);
	float sa = sin(Angle);

	tmp.Data[0][0] = 1; tmp.Data[0][1] = 0;  tmp.Data[0][2] = 0;   tmp.Data[0][3] = 0;
	tmp.Data[1][0] = 0; tmp.Data[1][1] = ca; tmp.Data[1][2] = -sa; tmp.Data[1][3] = 0;
	tmp.Data[2][0] = 0; tmp.Data[2][1] = sa; tmp.Data[2][2] = ca;  tmp.Data[2][3] = 0;
	tmp.Data[3][0] = 0; tmp.Data[3][1] = 0;  tmp.Data[3][2] = 0;   tmp.Data[3][3] = 1;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Set to a derivative matrix to rotate around the x axis.
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::GetdRotX(float Angle, float sign)
{
	Matrix3D tmp;

	float ca = cos(Angle);
	float sa = sin(Angle);

	tmp.Data[0][0] = 0; tmp.Data[0][1] = 0;  tmp.Data[0][2] = 0;   tmp.Data[0][3] = 0;
	tmp.Data[1][0] = 0; tmp.Data[1][1] = -sign*sa; tmp.Data[1][2] = -sign*ca; tmp.Data[1][3] = 0;
	tmp.Data[2][0] = 0; tmp.Data[2][1] = sign*ca; tmp.Data[2][2] = -sign*sa;  tmp.Data[2][3] = 0;
	tmp.Data[3][0] = 0; tmp.Data[3][1] = 0;  tmp.Data[3][2] = 0;   tmp.Data[3][3] = 0;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Set to a matrix to rotate around the y axis.
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::GetRotY(float Angle)
{
Matrix3D tmp;

	float ca = cos(Angle);
	float sa = sin(Angle);

	tmp.Data[0][0] = ca;  tmp.Data[0][1] = 0; tmp.Data[0][2] = sa; tmp.Data[0][3] = 0;
	tmp.Data[1][0] = 0;   tmp.Data[1][1] = 1; tmp.Data[1][2] = 0;  tmp.Data[1][3] = 0;
	tmp.Data[2][0] = -sa; tmp.Data[2][1] = 0; tmp.Data[2][2] = ca; tmp.Data[2][3] = 0;
	tmp.Data[3][0] = 0;   tmp.Data[3][1] = 0; tmp.Data[3][2] = 0;  tmp.Data[3][3] = 1;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Set to a derivative matrix to rotate around the y axis.
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::GetdRotY(float Angle, float sign)
{
	Matrix3D tmp;

	float ca = cos(Angle);
	float sa = sin(Angle);

	tmp.Data[0][0] = -sign*sa;  tmp.Data[0][1] = 0; tmp.Data[0][2] = sign*ca; tmp.Data[0][3] = 0;
	tmp.Data[1][0] = 0;   tmp.Data[1][1] = 0; tmp.Data[1][2] = 0;  tmp.Data[1][3] = 0;
	tmp.Data[2][0] = -sign*ca; tmp.Data[2][1] = 0; tmp.Data[2][2] = -sign*sa; tmp.Data[2][3] = 0;
	tmp.Data[3][0] = 0;   tmp.Data[3][1] = 0; tmp.Data[3][2] = 0;  tmp.Data[3][3] = 0;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////////////////
// Set to a matrix to rotate around the z axis.
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::GetRotZ(float Angle)
{
	Matrix3D tmp;

	float ca = cos(Angle);
	float sa = sin(Angle);

	tmp.Data[0][0] = ca; tmp.Data[0][1] = -sa; tmp.Data[0][2] = 0; tmp.Data[0][3] = 0;
	tmp.Data[1][0] = sa; tmp.Data[1][1] = ca;  tmp.Data[1][2] = 0; tmp.Data[1][3] = 0;
	tmp.Data[2][0] = 0;  tmp.Data[2][1] = 0;   tmp.Data[2][2] = 1; tmp.Data[2][3] = 0;
	tmp.Data[3][0] = 0;  tmp.Data[3][1] = 0;   tmp.Data[3][2] = 0; tmp.Data[3][3] = 1;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////////////////
// Set to a derivative matrix to rotate around the z axis.
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::GetdRotZ(float Angle, float sign)
{
	Matrix3D tmp;

	float ca = cos(Angle);
	float sa = sin(Angle);

	tmp.Data[0][0] = -sign*sa; tmp.Data[0][1] = -sign*ca; tmp.Data[0][2] = 0; tmp.Data[0][3] = 0;
	tmp.Data[1][0] = sign*ca; tmp.Data[1][1] = -sign*sa;  tmp.Data[1][2] = 0; tmp.Data[1][3] = 0;
	tmp.Data[2][0] = 0;  tmp.Data[2][1] = 0;   tmp.Data[2][2] = 0; tmp.Data[2][3] = 0;
	tmp.Data[3][0] = 0;  tmp.Data[3][1] = 0;   tmp.Data[3][2] = 0; tmp.Data[3][3] = 0;

	return tmp;
}


///////////////////////////////////////////////////////////////////////////
// Set to a matrix to perform a perspective transform with a focal spot
// at Dz on the z axis.
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::PerspectZ(float Dz)
{
Matrix3D tmp;

	tmp.Data[0][0] = 1; tmp.Data[0][1] = 0; tmp.Data[0][2] = 0;       tmp.Data[0][3] = 0;
	tmp.Data[1][0] = 0; tmp.Data[1][1] = 1; tmp.Data[1][2] = 0;       tmp.Data[1][3] = 0;
	tmp.Data[2][0] = 0; tmp.Data[2][1] = 0; tmp.Data[2][2] = 0;       tmp.Data[2][3] = 0;
	tmp.Data[3][0] = 0; tmp.Data[3][1] = 0; tmp.Data[3][2] = -1.0/Dz; tmp.Data[3][3] = 1;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Set to a matrix to perform a perspective transform with a focal spot
// at (Dx, Dy, Dz).
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::PerspectZ(float Dx, float Dy, float Dz)
{
	Matrix3D tmp;

	tmp.Data[0][0] = 1; tmp.Data[0][1] = 0; tmp.Data[0][2] = -Dx / Dz;       tmp.Data[0][3] = 0;
	tmp.Data[1][0] = 0; tmp.Data[1][1] = 1; tmp.Data[1][2] = -Dy / Dz;       tmp.Data[1][3] = 0;
	tmp.Data[2][0] = 0; tmp.Data[2][1] = 0; tmp.Data[2][2] = 0;              tmp.Data[2][3] = 0;
	tmp.Data[3][0] = 0; tmp.Data[3][1] = 0; tmp.Data[3][2] = -1.0/Dz;        tmp.Data[3][3] = 1;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Set to a matrix to perform a perspective transform with a focal spot
// at Dx on the x axis. Project onto YZ plane
///////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::PerspectX(float Dx)
{
    Matrix3D tmp;

	tmp.Data[0][0] = 0;       tmp.Data[0][1] = 0; tmp.Data[0][2] = 0;       tmp.Data[0][3] = 0;
	tmp.Data[1][0] = 0;       tmp.Data[1][1] = 1; tmp.Data[1][2] = 0;       tmp.Data[1][3] = 0;
	tmp.Data[2][0] = 0;       tmp.Data[2][1] = 0; tmp.Data[2][2] = 1;       tmp.Data[2][3] = 0;
	tmp.Data[3][0] = -1.0/Dx; tmp.Data[3][1] = 0; tmp.Data[3][2] = 0;       tmp.Data[3][3] = 1;

	tmp.Data[0][0] = 0;       tmp.Data[0][1] = 0; tmp.Data[0][2] = 0;       tmp.Data[0][3] = 0;
	tmp.Data[1][0] = 0;       tmp.Data[1][1] = -Dx; tmp.Data[1][2] = 0;       tmp.Data[1][3] = 0;
	tmp.Data[2][0] = 0;       tmp.Data[2][1] = 0; tmp.Data[2][2] = -Dx;       tmp.Data[2][3] = 0;
	tmp.Data[3][0] = 1.0;     tmp.Data[3][1] = 0; tmp.Data[3][2] = 0;       tmp.Data[3][3] = -Dx;

	return tmp;
}
///////////////////////////////////////////////////////////////////////////////////////////
// Set to a matrix to perform a perspective transform from a focal spot located at v = (xs,ys,zs,1)
// onto a plane described by n = (a,b,c,d) the plane vector of a plane described by the
// general equation ax + by + cz + d = 0.
// M = vn^T - <n,v>I I = identity matrix
/////////////////////////////////////////////////////////////////////////////////////////

Matrix3D Matrix3D::Perspect(float v[4], float n[4])
{
    Matrix3D tmp;

	float dot = v[0]*n[0] + v[1]*n[1] + v[2]*n[2] + v[3]*n[3];

	tmp.Data[0][0] = v[0]*n[0]; tmp.Data[0][1] = v[0]*n[1]; tmp.Data[0][2] = v[0]*n[2]; tmp.Data[0][3] = v[0]*n[3];
	tmp.Data[1][0] = v[1]*n[0]; tmp.Data[1][1] = v[1]*n[1]; tmp.Data[1][2] = v[1]*n[2]; tmp.Data[1][3] = v[1]*n[3];
	tmp.Data[2][0] = v[2]*n[0]; tmp.Data[2][1] = v[2]*n[1]; tmp.Data[2][2] = v[2]*n[2]; tmp.Data[2][3] = v[2]*n[3];
	tmp.Data[3][0] = v[3]*n[0]; tmp.Data[3][1] = v[3]*n[1]; tmp.Data[3][2] = v[3]*n[2]; tmp.Data[3][3] = v[3]*n[3];

	tmp.Data[0][0] -= dot;
	tmp.Data[1][1] -= dot;
	tmp.Data[2][2] -= dot;
	tmp.Data[3][3] -= dot;

	return tmp;
}

///////////////////////////////////////////////////////////////////////////
// Multiply Matrix with a vector (Dst = this * Src).
// Note:  For speed Src and Dst must point to different storage
///////////////////////////////////////////////////////////////////////////

void Matrix3D::VectorMult(float Src[4], float Dst[4]) const
{
	Dst[0] = Data[0][0]*Src[0] + Data[0][1]*Src[1] + Data[0][2]*Src[2] + Data[0][3]*Src[3];
	Dst[1] = Data[1][0]*Src[0] + Data[1][1]*Src[1] + Data[1][2]*Src[2] + Data[1][3]*Src[3];
	Dst[2] = Data[2][0]*Src[0] + Data[2][1]*Src[1] + Data[2][2]*Src[2] + Data[2][3]*Src[3];
	Dst[3] = Data[3][0]*Src[0] + Data[3][1]*Src[1] + Data[3][2]*Src[2] + Data[3][3]*Src[3];
}

///////////////////////////////////////////////////////////////////////////
// Transform a vector into homogeneous coordinates.
///////////////////////////////////////////////////////////////////////////

void Matrix3D::Homogenize(float Src[4], float Dst[4])
{
	Dst[0] = Src[0] / Src[3];
	Dst[1] = Src[1] / Src[3];
	Dst[2] = Src[2] / Src[3];
	Dst[3] = 1.0;
}

///////////////////////////////////////////////////////////////////////////
// Fast 2D projection (vector multiplication and homogenize).
// Note:  Makes assumptions about nature of matrix transform
///////////////////////////////////////////////////////////////////////////

void Matrix3D::Proj2D(float Src[3], float Dst[2]) const
{
register float h;

	h = Data[3][0]*Src[0] + Data[3][1]*Src[1] + Data[3][2]*Src[2] + Data[3][3];
	Dst[0] = (Data[0][0]*Src[0] + Data[0][1]*Src[1] + Data[0][2]*Src[2] + Data[0][3]) / h;
	Dst[1] = (Data[1][0]*Src[0] + Data[1][1]*Src[1] + Data[1][2]*Src[2] + Data[1][3]) / h;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//Convert rotation matrix to Euler angles : Rotation = GetRotZ() * GetRotY() * GetRotX()
////////////////////////////////////////////////////////////////////////////////////////////////////
void Matrix3D::RotationMatrixToEuler(float &angleX, float &angleY, float &angleZ)
{
	const float PI = 3.1415926536f;
	if(fabs(Data[2][0]) != 1.0f)
	{
		float angleY1 = asin(-Data[2][0]);
		float angleY2 = PI - angleY1;
		float angleX1 = atan2( Data[2][1] / cos(angleY1), Data[2][2] / cos(angleY1));
		float angleX2 = atan2( Data[2][1] / cos(angleY2), Data[2][2] / cos(angleY2));
		float angleZ1 = atan2( Data[1][0] / cos(angleY1), Data[0][0] / cos(angleY1));
		float angleZ2 = atan2( Data[1][0] / cos(angleY2), Data[0][0] / cos(angleY2));

		if( (angleX1*angleX1) < (angleX2*angleX2))
		{
			angleY = angleY1;
			angleX = angleX1;
			angleZ = angleZ1;
		}
		else
		{
			angleY = angleY2;
			angleX = angleX2;
			angleZ = angleZ2;
		}
	}
	else
	{
		angleZ = 0.0f;
		if(Data[2][0] == -1.0f)
		{
			angleY = PI / 2.0f;
			angleX = angleZ + atan2(Data[0][1], Data[0][2]);
		}
		else
		if(Data[2][0] == -1.0f)
		{
			angleY = -PI / 2.0f;
			angleX = -angleZ + atan2(Data[0][1], Data[0][2]);
		}
	}

	while(angleY < 0.0f) angleY += 2.0f * PI;
	while(angleX < -1.5f) angleX += 2.0f * PI;
	while(angleZ < -1.5f) angleZ += 2.0f * PI;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//return the matrix transpose
////////////////////////////////////////////////////////////////////////////////////////////////////
Matrix3D Matrix3D::GetTranspose(Matrix3D& Src)
{
	Matrix3D tmp;

	tmp.Data[0][0] = Src.Data[0][0]; tmp.Data[0][1] = Src.Data[1][0]; tmp.Data[0][2] = Src.Data[2][0]; tmp.Data[0][3] = Src.Data[3][0];
	tmp.Data[1][0] = Src.Data[0][1]; tmp.Data[1][1] = Src.Data[1][1]; tmp.Data[1][2] = Src.Data[2][1]; tmp.Data[1][3] = Src.Data[3][1];
	tmp.Data[2][0] = Src.Data[0][2]; tmp.Data[2][1] = Src.Data[1][2]; tmp.Data[2][2] = Src.Data[2][2]; tmp.Data[2][3] = Src.Data[3][2];
	tmp.Data[3][0] = Src.Data[0][3]; tmp.Data[3][1] = Src.Data[1][3]; tmp.Data[3][2] = Src.Data[2][3]; tmp.Data[3][3] = Src.Data[3][3];

	return tmp;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//return the inverse of [R|T]   matrix  = [RInv|RInv*T]
////////////////////////////////////////////////////////////////////////////////////////////////////
Matrix3D Matrix3D::GetInverseRT(Matrix3D& Src)
{
	Matrix3D tmp;

	//transpose = inverse rotational matrix
	tmp.Data[0][0] = Src.Data[0][0]; tmp.Data[0][1] = Src.Data[1][0]; tmp.Data[0][2] = Src.Data[2][0]; tmp.Data[0][3] = 0.0f;
	tmp.Data[1][0] = Src.Data[0][1]; tmp.Data[1][1] = Src.Data[1][1]; tmp.Data[1][2] = Src.Data[2][1]; tmp.Data[1][3] = 0.0f;
	tmp.Data[2][0] = Src.Data[0][2]; tmp.Data[2][1] = Src.Data[1][2]; tmp.Data[2][2] = Src.Data[2][2]; tmp.Data[2][3] = 0.0f;
	tmp.Data[3][0] = 0.0f;           tmp.Data[3][1] = 0.0f;           tmp.Data[3][2] = 0.0f;           tmp.Data[3][3] = 1.0f;

	float t[4], tI[4];
	t[0] = Src.Data[0][3];
	t[1] = Src.Data[1][3];
	t[2] = Src.Data[2][3];
	t[3] = Src.Data[3][3];

	tmp.VectorMult(t, tI);

	tmp.Data[0][3] = -tI[0];
	tmp.Data[1][3] = -tI[1];
	tmp.Data[2][3] = -tI[2];
	tmp.Data[3][3] =  1.0f;

	return tmp;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Invers of a 4x4 matrix (close form solution)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Matrix3D Matrix3D::GetInverse(Matrix3D& Src)
{
	Matrix3D tmp;
	float a11 = Src.Data[0][0];
	float a12 = Src.Data[0][1];
	float a13 = Src.Data[0][2];
	float a14 = Src.Data[0][3];

	float a21 = Src.Data[1][0];
	float a22 = Src.Data[1][1];
	float a23 = Src.Data[1][2];
	float a24 = Src.Data[1][3];

	float a31 = Src.Data[2][0];
	float a32 = Src.Data[2][1];
	float a33 = Src.Data[2][2];
	float a34 = Src.Data[2][3];

	float a41 = Src.Data[3][0];
	float a42 = Src.Data[3][1];
	float a43 = Src.Data[3][2];
	float a44 = Src.Data[3][3];

	float det = a11 * a22 * a33 * a44 + a11 * a23 * a34 * a42 + a11 * a24 * a32 * a43 - a11 * a24 * a33 * a42 - a11 * a23 * a32 * a44 - a11 * a22 * a34 * a43
		- a12 * a21 * a33 * a44 - a13 * a21 * a34 * a42 - a14 * a21 * a32 * a43 + a14 * a21 * a33 * a42 + a13 * a21 * a32 * a44 + a12 * a21 * a34 * a43
		+ a12 * a23 * a31 * a44 + a13 * a24 * a31 * a42 + a14 * a22 * a31 * a43 - a14 * a23 * a31 * a42 - a13 * a22 * a31 * a44 - a12 * a24 * a31 * a43
		- a12 * a23 * a34 * a41 - a13 * a24 * a32 * a41 - a14 * a22 * a33 * a41 + a14 * a23 * a32 * a41 + a13 * a22 * a34 * a41 + a12 * a24 * a33 * a41;

	tmp.Data[0][0] =  a22 * a33 * a44 + a23 * a34 * a42 + a24 * a32 * a43 - a24 * a33 * a42 - a23 * a32 * a44 - a22 * a34 * a43;
	tmp.Data[0][1] = -a12 * a33 * a44 - a13 * a34 * a42 - a14 * a32 * a43 + a14 * a33 * a42 + a13 * a32 * a44 + a12 * a34 * a43;
	tmp.Data[0][2] =  a12 * a23 * a44 + a13 * a24 * a42 + a14 * a22 * a43 - a14 * a23 * a42 - a13 * a22 * a44 - a12 * a24 * a43;
	tmp.Data[0][3] = -a12 * a23 * a34 - a13 * a24 * a32 - a14 * a22 * a33 + a14 * a23 * a32 + a13 * a22 * a34 + a12 * a24 * a33;

	tmp.Data[1][0] = -a21 * a33 * a44 - a23 * a34 * a41 - a24 * a31 * a43 + a24 * a33 * a41 + a23 * a31 * a44 + a21 * a34 * a43;
	tmp.Data[1][1] =  a11 * a33 * a44 + a13 * a34 * a41 + a14 * a31 * a43 - a14 * a33 * a41 - a13 * a31 * a44 - a11 * a34 * a43;
	tmp.Data[1][2] = -a11 * a23 * a44 - a13 * a24 * a41 - a14 * a21 * a43 + a14 * a23 * a41 + a13 * a21 * a44 + a11 * a24 * a43;
	tmp.Data[1][3] =  a11 * a23 * a34 + a13 * a24 * a31 + a14 * a21 * a33 - a14 * a23 * a31 - a13 * a21 * a34 - a11 * a24 * a33;

	tmp.Data[2][0] =  a21 * a32 * a44 + a22 * a34 * a41 + a24 * a31 * a42 - a24 * a32 * a41 - a22 * a31 * a44 - a21 * a34 * a42;
	tmp.Data[2][1] = -a11 * a32 * a44 - a12 * a34 * a41 - a14 * a31 * a42 + a14 * a22 * a41 + a12 * a31 * a44 + a11 * a34 * a42;
	tmp.Data[2][2] =  a11 * a22 * a44 + a12 * a24 * a41 + a14 * a21 * a42 - a14 * a22 * a41 - a12 * a21 * a44 - a11 * a24 * a42;
	tmp.Data[2][3] = -a11 * a22 * a34 - a12 * a24 * a31 - a14 * a21 * a32 + a14 * a22 * a31 + a12 * a21 * a34 + a11 * a24 * a32;

	tmp.Data[3][0] = -a21 * a32 * a43 - a22 * a33 * a41 - a23 * a31 * a42 + a23 * a32 * a41 + a22 * a31 * a43 + a21 * a33 * a42;
	tmp.Data[3][1] =  a11 * a32 * a43 + a12 * a33 * a41 + a13 * a31 * a42 - a13 * a32 * a41 - a12 * a31 * a43 - a11 * a33 * a42;
	tmp.Data[3][2] = -a11 * a22 * a43 - a12 * a23 * a41 - a13 * a21 * a42 + a13 * a22 * a41 + a12 * a21 * a43 + a11 * a23 * a42;
	tmp.Data[3][3] =  a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a13 * a22 * a31 - a12 * a21 * a33 - a11 * a23 * a32;

	float val = 1.0f / det;

	for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
			tmp.Data[j][i] *= val;


	return tmp;
}