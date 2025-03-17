#pragma once

struct Element
{
	int nodeIDs[4];
};

struct Mesh
{
	float* nodesX;
	float* nodesY;
	float* nodesZ;
	Element* elements;
	int nNodes;
	int nElements;
};

struct Constraint 
{
	int node;
	int DOF;
	float value;
};



struct extrinsicParms
{
	float tx;      // x translation  
	float ty;      // y translation
	float tz;      // z translation
	float thetaX;   //rotation about x axis
	float thetaY;     //rotation about y axis
	float thetaZ;     //rotation about z axis 
};

struct intrinsicParms
{
	float k0;      // lens parameter  
	float k1;      //lens parameter
	float p0;      // lens parameter  
	float p1;      //lens parameter
	float fx;      // focal length x
	float gamma;   //skew
	float uc;     // u center of sensro
	float fy;     //focal length y
	float vc;     //v center of sensor
};



struct Pt3D
{
	float x;
	float y;
	float z;
	float w;
};

struct Pt2D
{
	float x;
	float y;
	float w;
};

struct Point3D
{
	int X;
	int Y;
	int Z;
};

struct Point2D
{
	int X;
	int Y;
};

struct Geometry
{
	int ImageRows;
	int ImageCols;
	int VolumeSlices;
	int VolumeRows;
	int VolumeCols;
	float PixelSpacingX_cm;
	float PixelSpacingY_cm;
	float VoxelSpacingX_cm;
	float VoxelSpacingY_cm;
	float VoxelSpacingZ_cm;
};
