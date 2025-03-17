#include "pts_tools.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846  
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////
//  
//////////////////////////////////////////////////////////////////////////////////////////////////////
void rotateLiver()
{

	FILE* fiIn, *fiOut;
	int junk;

	fopen_s(&fiIn, "D:\\LiverMeshesDatabase\\Kelly_dataset\\FourLivers\\DLFlat_LM_732x732x443_Second.off", "r");
	fopen_s(&fiOut, "D:\\LiverMeshesDatabase\\Kelly_dataset\\FourLivers\\Flat_Surface.off", "w");

	//fscanf_s(fiIn, "33\n");
	fscanf_s(fiIn, "OFF\n");
	int nNodes, nElements;
	fscanf_s(fiIn, "%d %d %d\n", &nNodes, &nElements, &junk);

	fprintf_s(fiOut, "OFF\n");
	fprintf_s(fiOut, "%d %d %d\n", nNodes, nElements, 0);

	float nodeX, tempX;
	float nodeY, tempY;
	float nodeZ, tempZ;

	for (int i = 0; i < nNodes; i++)
	{
		//fscanf_s(fiIn, "%d %f %f %f\n", &junk, &nodeX, &nodeY, &nodeZ);
		fscanf_s(fiIn, "%f %f %f\n", &nodeX, &nodeY, &nodeZ);

		nodeX = abs(nodeX);
		nodeY = abs(nodeY);
		nodeZ = abs(nodeZ);
		tempX = 0.993113f * nodeX - 0.115917f * nodeY - 0.0170114f * nodeZ + 19.7226f;
		tempY = 0.115039f * nodeX + 0.992306f * nodeY - 0.04577372f * nodeZ - 12.3876f;
		tempZ = 0.0221865f * nodeX + 0.0435015f * nodeY + 0.998807f * nodeZ - 6.31452f;

		tempX = nodeX;
		tempY = nodeY;
		tempZ = nodeZ;

		fprintf_s(fiOut, "%f %f %f\n", tempX, tempY, tempZ);
	}

	int vert1, vert2, vert3, vert4;
	int fac1, fac2, fac3;
	for (int i = 0; i < nElements; i++)
	{
		//fscanf_s(fiIn, "%d %d %d %d %d %d %d\n", &fac1, &fac2, &fac2, &vert1, &vert2, &vert3, &vert4);
		//fprintf_s(fiOut, "%d %d %d %d %d\n", 4, vert1-1, vert2-1, vert3-1, vert4-1);
		fscanf_s(fiIn, "%d %d %d %d\n", &fac1, &vert1, &vert2, &vert3);
		fprintf_s(fiOut, "%d %d %d %d \n", 3, vert1, vert2, vert3);
	}

	fclose(fiIn);
	fclose(fiOut);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//  generate a partial surface of the 
//////////////////////////////////////////////////////////////////////////////////////////////////////
void generateIntraopSurface()
{

	FILE* fiIn, * fiOut;
	int junk;

	fopen_s(&fiIn, "D:\\LiverMeshesDatabase\\Simulations\\KellysPhantomLiver\\Sim1\\deformed_liver_Surface.off", "r");

	fscanf_s(fiIn, "OFF\n");
	int nNodes, nElements;
	fscanf_s(fiIn, "%d %d %d\n", &nNodes, &nElements, &junk);

	float* nodesX = new float[nNodes];
	float* nodesY = new float[nNodes];
	float* nodesZ = new float[nNodes];
	Element* elements = new Element[nElements];

	float node1;
	float node2;
	float node3;

	for (int i = 0; i < nNodes; i++)
	{
		fscanf_s(fiIn, "%f %f %f\n", &node1, &node2, &node3);

		nodesX[i] = node1;
		nodesY[i] = node2;
		nodesZ[i] = node3;
	}

	int vert1, vert2, vert3, vert4;
	int fac1, fac2, fac3;
	for (int i = 0; i < nElements; i++)
	{
		fscanf_s(fiIn, "%d %d %d %d\n", &fac1, &vert1, &vert2, &vert3);
		elements[i].nodeIDs[0] = vert1;
		elements[i].nodeIDs[1] = vert2;
		elements[i].nodeIDs[2] = vert3;
	}

	fclose(fiIn);

	Element* elementsOut = new Element[nElements];

	int cnt = 0;
	for (int i = 0; i < nElements; i++)
	{
		float Ax = nodesX[elements[i].nodeIDs[1]] - nodesX[elements[i].nodeIDs[0]];
		float Ay = nodesY[elements[i].nodeIDs[1]] - nodesY[elements[i].nodeIDs[0]];
		float Az = nodesZ[elements[i].nodeIDs[1]] - nodesZ[elements[i].nodeIDs[0]];
		float Bx = nodesX[elements[i].nodeIDs[2]] - nodesX[elements[i].nodeIDs[0]];
		float By = nodesY[elements[i].nodeIDs[2]] - nodesY[elements[i].nodeIDs[0]];
		float Bz = nodesZ[elements[i].nodeIDs[2]] - nodesZ[elements[i].nodeIDs[0]];

		float nz = Ax * By - Ay * Bx;
		float ny = Az * Bx - Ax * Bz;
		float nx = Ay * Bz - Az * By;

		float mag = sqrt(nx * nx + ny * ny + nz * nz);

		float ang = 180.0f * acos(nz/mag) / M_PI;
		if (ang < 60.0f)
		{
			elementsOut[cnt].nodeIDs[0] = elements[i].nodeIDs[0];
			elementsOut[cnt].nodeIDs[1] = elements[i].nodeIDs[1];
			elementsOut[cnt].nodeIDs[2] = elements[i].nodeIDs[2];
			cnt++;
		}
	}

	fopen_s(&fiOut, "D:\\LiverMeshesDatabase\\Simulations\\KellysPhantomLiver\\Sim1\\deformed_liver_intraop_surface_60.off", "w");
	fprintf_s(fiOut, "OFF\n");
	fprintf_s(fiOut, "%d %d %d\n", nNodes,cnt, 0);

	for (int i = 0; i < nNodes; i++)
	{
		fprintf_s(fiOut, "%f %f %f\n", nodesX[i], nodesY[i], nodesZ[i]);
	}

	for (int i = 0; i < cnt; i++)
	{
		fprintf_s(fiOut, "%d %d %d %d \n", 3, elementsOut[i].nodeIDs[0], elementsOut[i].nodeIDs[1], elementsOut[i].nodeIDs[2]);
	}

	fclose(fiOut);

	delete[] elements;
	delete[] elementsOut;
	delete[] nodesX;
	delete[] nodesY;
	delete[] nodesZ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//  
///////////////////////////////////////////////////////////////////////////////////////////////
void addGaussianNoise(Pt3D* ptCloud, int nPts, float sigmaX, float sigmaY, float sigmaZ)
{
	std::default_random_engine* generator = new std::default_random_engine(static_cast<long unsigned int>(time(0)));
	//std::default_random_engine* generator = new std::default_random_engine(static_cast<long unsigned int>(41151));
	std::normal_distribution<float>* distributionX = new std::normal_distribution<float>(0.0f, sigmaX);
	std::normal_distribution<float>* distributionY = new std::normal_distribution<float>(0.0f, sigmaY);
	std::normal_distribution<float>* distributionZ = new std::normal_distribution<float>(0.0f, sigmaZ);

	float del;
	float avg = 0.0f;
	float sd = 0.0f;
	for (int i = 0; i < nPts; i++)
	{
		ptCloud[i].x += (*distributionX)(*generator);
		ptCloud[i].y += (*distributionY)(*generator);
		ptCloud[i].z += (*distributionZ)(*generator);
		//del = (*distribution)(*generator);
		//avg += del;
		//sd += del * del;
	}

	//avg /= (float)nPts;
	//sd /= (float)nPts;

	//sd = sqrt(sd - avg * avg);

	delete generator;
	delete distributionX;
	delete distributionY;
	delete distributionZ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//  
///////////////////////////////////////////////////////////////////////////////////////////////
void rigidMotionTransformation(Pt3D * ptCloud, int nPts, extrinsicParms parms)
{

	float ptIn[4], ptOut[4];
	Matrix3D mat = Matrix3D::GetTrans(parms.tx, parms.ty, parms.tz) * Matrix3D::GetRotZ(parms.thetaZ) *
		                Matrix3D::GetRotY(parms.thetaY) * Matrix3D::GetRotX(parms.thetaX);
	for (int i = 0; i < nPts; i++)
	{
		ptIn[0] = ptCloud[i].x;
		ptIn[1] = ptCloud[i].y;
		ptIn[2] = ptCloud[i].z;
		ptIn[3] = 1.0f;

		mat.VectorMult(ptIn, ptOut);

		ptCloud[i].x = ptOut[0];
		ptCloud[i].y = ptOut[1];
		ptCloud[i].z = ptOut[2];
	}
}

float ScTP(const Pt3D& a, const Pt3D& b, const Pt3D& c)
{
	// computes scalar triple product

	//cross product b x c
	Pt3D d;
	d.x = b.y * c.z - b.z * c.y;
	d.y = b.z * c.x - b.x * c.z;
	d.z = b.x * c.y - b.y * c.x;

	//dot product a dot d
	float dot = a.x * d.x + a.y * d.y + a.z * d.z;

	return dot;
}

float VolTet(const Pt3D& a, const Pt3D& b, const Pt3D& c)
{
	// computes scalar triple product

	//cross product b x c
	Pt3D d;
	d.x = b.y * c.z - b.z * c.y;
	d.y = b.z * c.x - b.x * c.z;
	d.z = b.x * c.y - b.y * c.x;

	//dot product a dot d
	float dot = a.x * d.x + a.y * d.y + a.z * d.z;

	return 0.16667f * dot;
}

Pt3D crossProduct(const Pt3D& b, const Pt3D& c)
{
	//cross product b x c
	Pt3D d;
	d.x = b.y * c.z - b.z * c.y;
	d.y = b.z * c.x - b.x * c.z;
	d.z = b.x * c.y - b.y * c.x;

	return d;
}

Pt3D Subtract(const Pt3D& a, const Pt3D& b)
{
	Pt3D dif;
	dif.x = b.x - a.x;
	dif.y = b.y - a.y;
	dif.z = b.z - a.z;
	dif.w = b.w - a.w;
	return dif;
}

void ConvertXYZtoOFF()
{
	FILE* fiIn, *fiOut;
	string fileIn = "D:\\LiverMeshesDatabase\\Sparse_data_challenge\\AdditionalDataForParticipants\\Set084_TargetsIncomplete.xyz";
	string fileOut = "D:\\LiverMeshesDatabase\\Sparse_data_challenge\\AdditionalDataForParticipants\\Set084_TargetsIncomplete.off";

	int nPts = 35;

	fopen_s(&fiIn, fileIn.c_str(), "r");
	fopen_s(&fiOut, fileOut.c_str(), "w");

	fprintf_s(fiOut, "OFF\n");
	fprintf_s(fiOut, " %d, 0, 0\n", nPts);

	int junk;
	float x, y, z;
	for (int i = 0; i < nPts; i++)
	{
		fscanf_s(fiIn,"%d %f %f %f\n", &junk, &x, &y, &z);
		fprintf_s(fiOut, "%f %f %f\n", x, y, z);
	}

	fclose(fiOut);
	fclose(fiIn);
}

